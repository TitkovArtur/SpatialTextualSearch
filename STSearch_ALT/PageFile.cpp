#include "PageFile.h"

PageFile::PageFile(const char *name, int p_len, const char *path, c_policy pol, int c_size, bool force_new)
{
	EmptyValue();

	Construct(name, p_len, path, pol, c_size, force_new);
}
PageFile::PageFile(int p_len, PageFile::c_policy pol, int c_size)
{
    EmptyValue();

    Construct(p_len, pol, c_size);
}
PageFile::~PageFile()
{
	CleanValue();
}
void PageFile::Construct(int p_len, c_policy pol, int c_size)
{
    //assign values
    page_len=p_len;
    page_num=1; //we always have a header page
    policy=pol;
    if((pol==C_FULLMEM)||(pol==C_NO_CACHE)) cache_size=0;
    else cache_size=c_size;

    //we then allocate one buffer page to use
    buf_page=new char[page_len];    //page_len = 4096
    memset(buf_page, 0, page_len);
    memcpy(buf_page, &page_len, sizeof(int));
    memcpy(buf_page+sizeof(int), &page_num, sizeof(int));

    switch(policy)
    {
    case C_FULLMEM:
        cache=new vector<char *>();
        cache->resize(page_num, NULL);
        c_dirty=new vector<bool>();
        c_dirty->resize(page_num, false);

        //make everything cached
        for(int i=0;i<page_num;i++)    //page_num == 1
        {
            cache->at(i)=buf_page;
            c_dirty->at(i)=false;
        }

        break;
    default:
        printf("Wrong policy!\n");
        exit(1);
        break;
    }
}

void PageFile::Construct(const char *name, int p_len, const char *path, c_policy pol, int c_size, bool force_new)
{
	//assign values
	int fn_len=GenericTool::GetCombinedPath(path, name, NULL);
	fname=new char[fn_len];
	GenericTool::GetCombinedPath(path, name, fname);
	page_len=p_len;
	page_num=1; //we always have a header page
	policy=pol;
	if((pol==C_FULLMEM)||(pol==C_NO_CACHE)) cache_size=0;
	else cache_size=c_size;

	//test if we can find the file
	if(GenericTool::JudgeExistence(fname, force_new))
	{
		//file exists, then we read parameter from file
		f=new fstream(fname, ios::in | ios::out | ios::binary);
		new_file=false;

		if(!f->is_open())
		{
			printf("Open File %s Failed!\n", fname);
			exit(-1);
		}

		f->read((char *)&page_len, sizeof(int));
		f->read((char *)&page_num, sizeof(int));

		//allocate one buffer page to use
		buf_page=new char[page_len];
	}
	else
	{
		//file does not exist, we first create a zero-sized file by specify ios::out
		//then we close it and re-open it with ios::in and ios::out
		//because directly go on with only ios::out makes us cannot read
		f=new fstream(fname, ios::in | ios::out | ios::binary);
		new_file=true;

		if(!f->is_open())
		{
			printf("Open File %s Failed!\n", fname);
			exit(-1);
		}

		//we then allocate one buffer page to use, and write data to file
		buf_page=new char[page_len];
		memcpy(buf_page, &page_len, sizeof(int));
		memcpy(buf_page+sizeof(int), &page_num, sizeof(int));
		f->write(buf_page, page_len);
	}

	char *buf=NULL;

	if(f->fail()) printf("fail start!");

	switch(policy)
	{
	case C_FULLMEM:
		cache=new vector<char *>();
		cache->resize(page_num, NULL);
		c_dirty=new vector<bool>();
		c_dirty->resize(page_num, false);

		//make everything cached
		f->seekg(0, ios::beg);
		for(int i=0;i<page_num;i++)
		{
			buf=new char[page_len];
			f->read(buf, page_len);
			cache->at(i)=buf;
			c_dirty->at(i)=false;
		}

		break;
	case C_LRU:
	case C_MRU:
		cache=new vector<char *>();
		cache->resize(page_num, NULL);
		c_dirty=new vector<bool>();
		c_dirty->resize(page_num, false);
		use_list=new list<int>();
		list_pos=new vector<list<int>::iterator>();
		list_pos->resize(page_num, use_list->end());

		//make header cached
		buf=new char[page_len];
		f->seekg(0, ios::beg);
		f->read(buf, page_len);
		cache->at(0)=buf;
		c_dirty->at(0)=false;
		use_list->push_front(0);
		list_pos->at(0)=use_list->begin();

		cache_size=max(cache_size, 1);
		num_cached=1;

		break;
	}
}

void PageFile::EmptyValue()
{
	f=NULL;
	fname=NULL;
	new_file=false;
	page_len=0;
	page_num=0;
	policy=C_FULLMEM;
	cache_size=0;
	buf_page=NULL;
	num_cached=0;
	cache=NULL;
	c_dirty=NULL;
	use_list=NULL;
	list_pos=NULL;
	page_access=0;
}

void PageFile::CleanValue()
{
	flush();

	if(f!=NULL) f->close();
	if(fname!=NULL) delete[] fname;
	if(buf_page!=NULL) delete[] buf_page;
	if(cache!=NULL)
	{
		for(vector<char *>::iterator iter=cache->begin();iter!=cache->end();++iter)
		{
			delete[] (*iter);
		}
		delete cache;
	}
	if(c_dirty!=NULL) delete c_dirty;
	if(use_list!=NULL) delete use_list;
	if(list_pos!=NULL) delete list_pos;

	EmptyValue();
}

const char *PageFile::read_header()
{
	//return the address of page 0 plus the size of page_len+page_num
	return read_page(-1)+2*sizeof(int);
}

void PageFile::set_header(const char *header)
{
	memcpy(buf_page, &page_len, sizeof(int));
	memcpy(buf_page+sizeof(int), &page_num, sizeof(int));
	memcpy(buf_page+2*sizeof(int), header, page_len-2*sizeof(int));
	write_page(buf_page, -1);
}

const char *PageFile::read_page(int index)
{
	index++; //transfer to correct page index
	if((index<0)||(index>=page_num)) return NULL;

	//statistics
	page_access++;

	char *buf=NULL;
	int evict_ind=-1;

	switch(policy)
	{
	case C_FULLMEM:
		return cache->at(index);
		break;
	case C_LRU:
		if(cache->at(index)==NULL)
		{
			//no cache hit!
			if(num_cached<cache_size)
			{
				//cache underfull
				buf=new char[page_len];
				num_cached++;
			}
			else
			{
				//cache full, evict one page
				evict_ind=use_list->back();
				buf=cache->at(evict_ind);
				if(c_dirty->at(evict_ind))
				{
					//if dirty then we nned to write back to disk
					f->seekp(evict_ind*page_len, ios::beg);
					f->write(buf, page_len);
					c_dirty->at(evict_ind)=false;
				}
				cache->at(evict_ind)=NULL;
				use_list->pop_back();
				list_pos->at(evict_ind)=use_list->end();
			}

			//we cache new page
			f->seekg(index*page_len, ios::beg);
			f->read(buf, page_len);
			cache->at(index)=buf;
			c_dirty->at(index)=false;
			use_list->push_front(index);
			list_pos->at(index)=use_list->begin();
		}
		else
		{
			buf=cache->at(index);
			use_list->erase(list_pos->at(index));
			use_list->push_front(index);
			list_pos->at(index)=use_list->begin();
		}

		return buf;
		break;
	case C_MRU:
		if(cache->at(index)==NULL)
		{
			//no cache hit!
			if(num_cached<cache_size)
			{
				//cache underfull
				buf=new char[page_len];
				num_cached++;
			}
			else
			{
				//cache full, evict one page
				evict_ind=use_list->front();
				buf=cache->at(evict_ind);
				if(c_dirty->at(evict_ind))
				{
					//if dirty then we nned to write back to disk
					f->seekp(evict_ind*page_len, ios::beg);
					f->write(buf, page_len);
					c_dirty->at(evict_ind)=false;
				}
				cache->at(evict_ind)=NULL;
				use_list->pop_front();
				list_pos->at(evict_ind)=use_list->end();
			}

			//we cache new page
			f->seekg(index*page_len, ios::beg);
			f->read(buf, page_len);
			cache->at(index)=buf;
			c_dirty->at(index)=false;
			use_list->push_front(index);
			list_pos->at(index)=use_list->begin();
		}
		else
		{
			buf=cache->at(index);
			use_list->erase(list_pos->at(index));
			use_list->push_front(index);
			list_pos->at(index)=use_list->begin();
		}

		return buf;
		break;
	case C_NO_CACHE:
		f->seekg(index*page_len, ios::beg);
		f->read(buf_page, page_len);
		return buf_page;
		break;
	}

	return NULL;
}

bool PageFile::write_page(const char *p, int index)
{
	index++; //transfer to correct page index
	if((index<0)||(index>=page_num)) return false;

	//statistics
	page_access++;

	char *buf=NULL;
	int evict_ind=-1;

	switch(policy)
	{
	case C_FULLMEM:
		memcpy(cache->at(index), p, page_len);
		c_dirty->at(index)=true;
		break;
	case C_LRU:
		if(cache->at(index)==NULL)
		{
			//no cache hit!
			if(num_cached<cache_size)
			{
				//cache underfull
				buf=new char[page_len];
				num_cached++;
			}
			else
			{
				//cache full, evict one page
				evict_ind=use_list->back();
				buf=cache->at(evict_ind);
				if(c_dirty->at(evict_ind))
				{
					//if dirty then we nned to write back to disk
					f->seekp(evict_ind*page_len, ios::beg);
					f->write(buf, page_len);
					c_dirty->at(evict_ind)=false;
				}
				cache->at(evict_ind)=NULL;
				use_list->pop_back();
				list_pos->at(evict_ind)=use_list->end();
			}

			//we cache new page
			memcpy(buf, p, page_len);
			cache->at(index)=buf;
			c_dirty->at(index)=true;
			use_list->push_front(index);
			list_pos->at(index)=use_list->begin();
		}
		else
		{
			memcpy(cache->at(index), p, page_len);
			c_dirty->at(index)=true;
			use_list->erase(list_pos->at(index));
			use_list->push_front(index);
			list_pos->at(index)=use_list->begin();
		}
		break;
	case C_MRU:
		if(cache->at(index)==NULL)
		{
			//no cache hit!
			if(num_cached<cache_size)
			{
				//cache underfull
				buf=new char[page_len];
				num_cached++;
			}
			else
			{
				//cache full, evict one page
				evict_ind=use_list->front();
				buf=cache->at(evict_ind);
				if(c_dirty->at(evict_ind))
				{
					//if dirty then we nned to write back to disk
					f->seekp(evict_ind*page_len, ios::beg);
					f->write(buf, page_len);
					c_dirty->at(evict_ind)=false;
				}
				cache->at(evict_ind)=NULL;
				use_list->pop_front();
				list_pos->at(evict_ind)=use_list->end();
			}

			//we cache new page
			memcpy(buf, p, page_len);
			cache->at(index)=buf;
			c_dirty->at(index)=true;
			use_list->push_front(index);
			list_pos->at(index)=use_list->begin();
		}
		else
		{
			memcpy(cache->at(index), p, page_len);
			c_dirty->at(index)=true;
			use_list->erase(list_pos->at(index));
			use_list->push_front(index);
			list_pos->at(index)=use_list->begin();
		}
		break;
	case C_NO_CACHE:
		f->seekp(index*page_len, ios::beg);
		f->write(p, page_len);
		break;
	}

	return true;
}

int PageFile::append_page(const char *p)
{
	int index=page_num++;

	//statistics
	page_access++;

	char *buf=NULL;
	int evict_ind=-1;

	switch(policy)
	{
	case C_FULLMEM:
		buf=new char[page_len];
		memcpy(buf, p, page_len);
		cache->push_back(buf);
		c_dirty->push_back(true);
		break;
	case C_LRU:
		if(num_cached<cache_size)
		{
			//cache underfull
			buf=new char[page_len];
			num_cached++;
		}
		else
		{
			//cache full, evict one page
			evict_ind=use_list->back();
			buf=cache->at(evict_ind);
			if(c_dirty->at(evict_ind))
			{
				//if dirty then we nned to write back to disk
				f->seekp(evict_ind*page_len, ios::beg);
				f->write(buf, page_len);
				c_dirty->at(evict_ind)=false;
			}
			cache->at(evict_ind)=NULL;
			use_list->pop_back();
			list_pos->at(evict_ind)=use_list->end();
		}

		//we cache new page
		memcpy(buf, p, page_len);
		cache->push_back(buf);
		c_dirty->push_back(true);
		use_list->push_front(index);
		list_pos->push_back(use_list->begin());
		break;
	case C_MRU:
		if(num_cached<cache_size)
		{
			//cache underfull
			buf=new char[page_len];
			num_cached++;
		}
		else
		{
			//cache full, evict one page
			evict_ind=use_list->front();
			buf=cache->at(evict_ind);
			if(c_dirty->at(evict_ind))
			{
				//if dirty then we nned to write back to disk
				f->seekp(evict_ind*page_len, ios::beg);
				f->write(buf, page_len);
				c_dirty->at(evict_ind)=false;
			}
			cache->at(evict_ind)=NULL;
			use_list->pop_front();
			list_pos->at(evict_ind)=use_list->end();
		}

		//we cache new page
		memcpy(buf, p, page_len);
		cache->push_back(buf);
		c_dirty->push_back(true);
		use_list->push_front(index);
		list_pos->push_back(use_list->begin());
		break;
	case C_NO_CACHE:
		f->seekp(index*page_len, ios::beg);
		f->write(p, page_len);
		break;
	}

	return index-1;
}


//this is a cross-platform implementations
//since this function is slow, so avoid 
bool PageFile::truncate_pages(int num)
{
	//make sure all things are written back to disk
	flush();

	//release all memory cache
	if(cache!=NULL)
	{
		for(vector<char *>::iterator iter=cache->begin();iter!=cache->end();++iter)
		{
			delete[] (*iter);
		}
		delete cache;
	}
	if(c_dirty!=NULL) delete c_dirty;
	if(use_list!=NULL) delete use_list;
	if(list_pos!=NULL) delete list_pos;

	//reset values
	num_cached=0;
	page_num=max(page_num-num, 1);

	//open temporary file
	fstream *ftemp=new fstream("temp_page_file", ios::out | ios::binary);

	//copy file in pages
	f->seekg(0, ios::beg);
	ftemp->seekp(0, ios::beg);
	for(int i=0;i<page_num;i++)
	{
		f->read(buf_page, page_len);
		ftemp->write(buf_page, page_len);
	}

	//close both file
	f->close();
	delete f;
	ftemp->close();
	delete ftemp;

	//rename file
	if(remove(fname)||rename("temp_page_file", fname))
	{
		printf("Truncate file fail! Page file object is corrupted!\n");
		return false;
	}

	//reconstruct data structure
	char *buf=NULL;

	f=new fstream(fname, ios::in | ios::out | ios::binary);
	switch(policy)
	{
	case C_FULLMEM:
		cache=new vector<char *>();
		cache->resize(page_num, NULL);
		c_dirty=new vector<bool>();
		c_dirty->resize(page_num, false);

		//make everything cached
		f->seekg(0, ios::beg);
		for(int i=0;i<page_num;i++)
		{
			buf=new char[page_len];
			f->read(buf, page_len);
			cache->at(i)=buf;
			c_dirty->at(i)=false;
		}

		break;
	case C_LRU:
	case C_MRU:
		cache=new vector<char *>();
		cache->resize(page_num, NULL);
		c_dirty=new vector<bool>();
		c_dirty->resize(page_num, false);
		use_list=new list<int>();
		list_pos=new vector<list<int>::iterator>();
		list_pos->resize(page_num, use_list->end());

		//make header cached
		buf=new char[page_len];
		f->seekg(0, ios::beg);
		f->read(buf, page_len);
		cache->at(0)=buf;
		c_dirty->at(0)=false;
		use_list->push_front(0);
		list_pos->at(0)=use_list->begin();

		num_cached=1;

		break;
	}

	return true;
}

void PageFile::flush()
{
	switch(policy)
	{
	case C_FULLMEM:
		//check if we can flush
		if((c_dirty==NULL)||(cache==NULL)) return;

		for(int i=0;i<page_num;i++)
		{
			if(c_dirty->at(i))
			{
				c_dirty->at(i)=false;
				f->seekp(i*page_len, ios::beg);
				f->write(cache->at(i), page_len);
			}
		}
		break;
	case C_LRU:
	case C_MRU:
		//check if we can flush
		if((use_list==NULL)||(c_dirty==NULL)||(cache==NULL)) return;

		for(list<int>::iterator iter=use_list->begin();iter!=use_list->end();++iter)
		{
			int ind=*iter;
			if(c_dirty->at(ind))
			{
				c_dirty->at(ind)=false;
				f->seekp(ind*page_len, ios::beg);
				f->write(cache->at(ind), page_len);
			}
		}
		break;
	case C_NO_CACHE:
		break;
	}

	//make page_len and page_num synchronized
	f->seekp(0, ios::beg);
	f->write((char *)&page_len, sizeof(int));
	f->write((char *)&page_num, sizeof(int));
	f->flush();
}

/*bool PageFile::testPageFile(ConfigTool &cfg)
{
	PageFile::c_policy pols[]={PageFile::C_FULLMEM, PageFile::C_LRU, PageFile::C_MRU, PageFile::C_NO_CACHE};

	PageFile pf(cfg.GetConfigStr("filename"), 16, NULL, pols[cfg.GetConfigInt("cache_policy")],
		cfg.GetConfigInt("cache_size"), (cfg.GetConfigInt("force_new")!=0));

	char header[]="BUFFHEAD";
	char page[16];
	int twpage[5]={5,2,1,7,6};
	char final_content[]="\x10\0\0\0\x9\0\0\0BUFFHEAD0000000000000000xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx3333333333333333"
		"4444444444444444xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx";

	pf.set_header(header);
	for(int i=0;i<16;i++)
	{
		for(int j=0;j<16;j++) page[j]='0'+i;
		pf.append_page(page);
	}

	pf.truncate_pages(8);

	for(int j=0;j<16;j++) page[j]='x';

	for(int i=0;i<5;i++) pf.write_page(page, twpage[i]);

	pf.flush(); //means we destruct

	char buf[144];

	char *buffer=new char[GenericTool::GetCombinedPath(NULL, cfg.GetConfigStr("filename"), NULL)];
	GenericTool::GetCombinedPath(NULL, cfg.GetConfigStr("filename"), buffer);

	ifstream the_pf(buffer, ios::binary);

	if(!the_pf.is_open()) return false;

	the_pf.read(buf, 144);

	if(the_pf.fail()) return false;

	the_pf.close();

	for(int i=0;i<144;i++)
	{
		if(buf[i]!=final_content[i])
		{
			return false;
		}
	}

	delete[] buffer;

	return true;
}*/
