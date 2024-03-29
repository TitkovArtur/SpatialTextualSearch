#include "GenericTool.h"

//Cross Platform snprintf
#ifdef _MSC_VER
//Under vc, we have to use some simulation
int msvc_snprintf(char *str, size_t size, const char *format, ...)
{
	int count;

	va_list ap;
	va_start(ap, format);
	count=_vscprintf(format, ap);
	_vsnprintf_s(str, size, _TRUNCATE, format, ap);
	va_end(ap);

	return count;
}
#else
#ifdef __GNUC__
//We have already use the snprintf directly
#else
//For other compiler, we output error
int other_snprintf(char *str, size_t size, const char *format, ...)
{
	printf("Not Implemented!\n");

	return -1;
}
#endif
#endif

//file related tool functions
bool GenericTool::CheckPathExistence(const char *path)
{
	errno=0;
	return (!rename(path, path) || (errno!=2)); //can rename or errno is not "file not exist"
}

//if buffer is NULL, thren return the number of buffer required including '\0'
int GenericTool::RegularizeDirPath(const char *path, char *buffer)
{
	//default path assignment
#ifdef WIN32
	if(path==NULL) path="../"; //win32 default
#else
	if(path==NULL) path="./"; //linux default
#endif	

	int len=strlen(path);

	if(len==0) //path=""
	{
		if(buffer!=NULL) buffer[0]='\0';
		return 1;
	}

	//only check last character is '\\' or '/'
	if((path[len-1]=='\\')||(path[len-1]=='/'))
	{
		//OK, just copy
		if(buffer!=NULL)
		{
			memcpy(buffer, path, (len+1)*sizeof(char));
			
			//since UNIX does not recognize '\\' in path but windows do recognize '/', we replace all '\\' to '/'
			for(int i=0;i<len;i++) if(buffer[i]=='\\') buffer[i]='/';
		}			
	}
	else
	{
		//not OK, need one more
		len++;
		if(buffer!=NULL)
		{
			c99_snprintf(buffer, len+1, "%s/", path); //use linux default '/', generally windows can recognize

			//since UNIX does not recognize '\\' in path but windows do recognize '/', we replace all '\\' to '/'
			for(int i=0;i<len;i++) if(buffer[i]=='\\') buffer[i]='/';
		}
	}
	
	return len+1;
}

void GenericTool::EnsurePathExistence(const char *path)
{
	if((path!=NULL)&&!CheckPathExistence(path))
	{
		int len=strlen(path);
		char *dir_path=new char[len+10];
#ifdef WIN32
		int ret=c99_snprintf(dir_path, len+10, "mkdir \"%s\"", path); // win32
#else
		int ret=c99_snprintf(dir_path, len+10, "mkdir %s -p", path); // unix
#endif		
		system(dir_path); //use mkdir, which is suported by a variety of OS
		delete[] dir_path;
	}
}

//if buffer is NULL, thren return the number of buffer required including '\0'
int GenericTool::GetCombinedPath(const char *dir, const char *file, char *buffer)
{
	int len_dir=RegularizeDirPath(dir, buffer);
	int len_file=strlen(file);
	int len=len_dir+len_file;

	if(buffer!=NULL) memcpy(buffer+len_dir-1, file, (len_file+1)*sizeof(char));

	return len;
}

//if we are forced to create new file, then we will remove orginal file and return false
bool GenericTool::JudgeExistence(const char *full_path, bool force_new)
{
	bool new_buffer=false;
	int path_len=strlen(full_path);

	char *buffer=new char[path_len+1];
	memcpy(buffer, full_path, (path_len+1)*sizeof(char));
	for(int i=0;i<path_len;i++) if(buffer[i]=='\\') buffer[i]='/';

	if(CheckPathExistence(buffer)&&!force_new)
	{
		delete[] buffer;
		return true;
	}
	else
	{
		errno=0;
		if(remove(buffer)&&errno!=2) //if cannot remove and the reason is not that file not exist
		{
			printf("Cannot remove file: %s for a forced new file creation.\n", buffer);
			exit(-1);
		}
		if(force_new)
		{
			int ind;
			for(ind=path_len-1;ind>=0;ind--) if(buffer[ind]=='/') break;
			ind++;
			char temp=buffer[ind];
			buffer[ind]='\0';
			EnsurePathExistence(buffer);
			buffer[ind]=temp;
			fstream f(buffer, ios::out | ios::trunc); //make sure file exists for later open with zero size
			f.close();
		}
		delete[] buffer;
		return false;
	}
}

int GenericTool::ChangeFileExtension(const char *full_path, const char *new_ext, char *buffer)
{
	int o_len=strlen(full_path);
	int ext_len=strlen(new_ext);

	int pos=o_len-1;
	while((pos>=0)&&(full_path[pos]!='\\')&&(full_path[pos]!='/')&&(full_path[pos]!='.')) pos--;
	if(full_path[pos]!='.')
	{
		if(buffer!=NULL)
		{
			memcpy(buffer, full_path, o_len*sizeof(char));
			buffer[o_len]='.';
		}
		pos=o_len+1;
	}
	else
	{
		pos++;
		if(buffer!=NULL) memcpy(buffer, full_path, pos*sizeof(char));
	}

	if(buffer!=NULL) memcpy(buffer+pos, new_ext, (ext_len+1)*sizeof(char));

	return pos+ext_len+1;
}

//Some Cross Platform Important Functions
double GenericTool::GetGaussianRandom(double mean, double sigma)
{
	double v1, v2, s;

	do
	{
		v1=2*getrand()-1;
		v2=2*getrand()-1;
		s=v1*v1+v2*v2;
	} while(s>=1.0);

	if(s==0.0) return mean;
	else return v1*sqrt(-2.0*log(s)/s)*sigma+mean;
}

void GenericTool::Independent(float *data_f, int dim)
{
	for(int i=0;i<dim;i++) data_f[i]=(float)getrand();
}

void GenericTool::Independent(double *data_d, int dim)
{
	for(int i=0;i<dim;i++) data_d[i]=getrand();
}

void GenericTool::Correlated(float *data_f, int dim, float sigma)
{
	float value=(float)getrand();

	for(int i=0;i<dim;i++)
	{
		float newval;
		do
		{
			newval=(float)GetGaussianRandom(value, sigma);
		} while((newval<0.0f)||(newval>1.0f));
		data_f[i]=newval;
	}
}

void GenericTool::Correlated(double *data_d, int dim, double sigma)
{
	double value=getrand();

	for(int i=0;i<dim;i++)
	{
		double newval;
		do
		{
			newval=GetGaussianRandom(value, sigma);
		} while((newval<0.0)||(newval>1.0));
		data_d[i]=newval;
	}
}

//we generate every dim a random number in [0,1), and then normalize so that the sum equal to range
//then if any dim is out of range, we measure its distance with (range/dim, ..., range/dim) represented by the ratio
//we get the minimum ratio so that we can ensure all dim is in thr range
//in order to avoid lot of border data, we multiply the ratio by a random number in [0,1)
//so the final point is (range/dim, ..., range/dim)...ratio...final point...1-ratio...(normalized point)
void GenericTool::AntiCorrelated(float *data_f, int dim, float sigma)
{
	float range=max(min(0.5f*dim+(float)GetGaussianRandom(0, sigma), 1.0f*dim), 0.0f);

	float sum=0.0f, l_min=1.0f;
	for(int i=0;i<dim;i++) sum+=data_f[i]=(float)getrand();
	for(int i=0;i<dim;i++)
	{
		data_f[i]*=range/sum;
		if(data_f[i]>1.0f)
		{
			float l=(1.0f-range/dim)/(data_f[i]-range/dim);
			if(l<l_min) l_min=l;
		}
	}
	if(l_min<1.0f)
	{
		l_min*=(float)getrand();
		for(int i=0;i<dim;i++) data_f[i]=(1.0f-l_min)*range/dim+l_min*data_f[i];
	}
}

void GenericTool::AntiCorrelated(double *data_d, int dim, double sigma)
{
	double range=max(min(0.5*dim+GetGaussianRandom(0, sigma), 1.0*dim), 0.0);

	double sum=0.0, l_min=1.0;
	for(int i=0;i<dim;i++) sum+=data_d[i]=getrand();
	for(int i=0;i<dim;i++)
	{
		data_d[i]*=range/sum;
		if(data_d[i]>1.0)
		{
			double l=(1.0-range/dim)/(data_d[i]-range/dim);
			if(l<l_min) l_min=l;
		}
	}
	if(l_min<1.0)
	{
		l_min*=getrand();
		for(int i=0;i<dim;i++) data_d[i]=(1.0-l_min)*range/dim+l_min*data_d[i];
	}
}

void GenericTool::Clustered(float *data_f, int dim, float **centers, int num_center, float sigma)
{
	int cur_center=(int)(getrand()*num_center);

	for(int i=0;i<dim;i++)
	{
		float newval;
		do
		{
			newval=(float)GetGaussianRandom(centers[cur_center][i], sigma);
		} while((newval<0.0f)||(newval>1.0f));
		data_f[i]=newval;
	}
}

void GenericTool::Clustered(double *data_d, int dim, double **centers, int num_center, double sigma)
{
	int cur_center=(int)(getrand()*num_center);

	for(int i=0;i<dim;i++)
	{
		double newval;
		do
		{
			newval=GetGaussianRandom(centers[cur_center][i], sigma);
		} while((newval<0.0)||(newval>1.0));
		data_d[i]=newval;
	}
}