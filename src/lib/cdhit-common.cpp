// =============================================================================
// CD-HI/CD-HIT
//
// Cluster Database at High Identity
//
// CD-HIT clusters protein sequences at high identity threshold.
// This program can remove the high sequence redundance efficiently.
//
// program written by
//                    Weizhong Li
//                    UCSD, San Diego Supercomputer Center
//                    La Jolla, CA, 92093
//                    Email liwz@sdsc.edu
//
//                 at
//                    Adam Godzik's lab
//                    The Burnham Institute
//                    La Jolla, CA, 92037
//                    Email adam@burnham-inst.org
//
// modified by:
//                    Limin Fu
//                    Center for Research in Biological Systems (CRBS), UCSD
//                    La Jolla, CA, 92093
//                    Email: l2fu@ucsd.edu, fu@daovm.net
// =============================================================================

import "cdhit-common.h";

import Options;

#include<valarray>
#include<stdint.h>
#include<assert.h>
#include<limits.h>
#include<assert.h>

struct TempFile
{
	FILE *file;
	char buf[512];

	TempFile( const char *dir = NULL ){
		int len = dir ? strlen( dir ) : 0;
		assert( len < 400 );
		buf[0] = 0;
		if( len ){
			strcat( buf, dir );
			if( buf[len-1] != '/' && buf[len-1] != '\\' ){
				buf[len] = '/';
				len += 1;
			}
		}
		strcat( buf, "cdhit.temp." );
		len += 11;
		sprintf( buf + len, "%p", this );
		file = fopen( buf, "w+" );
	}
	~TempFile(){
		if( file ){
			fclose( file );
			remove( buf );
		}
	}
};

struct TempFiles
{
	NVector<TempFile*> files;

	~TempFiles(){ Clear(); }

	void Clear(){
		int i;
		#pragma omp critical
		{
			for(i=0; i<files.size; i++) if( files[i] ) delete files[i];
			files.Clear();
		}
	}
};

TempFiles temp_files;

FILE* OpenTempFile( const char *dir = NULL )
{
	TempFile *file = new TempFile( dir );
	#pragma omp critical
	{
		temp_files.files.push_back( file );
	}
	return file->file;
}

std::vector<size_t> Comp_AAN_idx;

void make_comp_iseq(int len, char *iseq_comp, char *iseq)
{
	int i, c[6] = {3,2,1,0,4,5};
	for (i=0; i<len; i++) iseq_comp[i] = c[ (int)iseq[len-i-1] ];
}

void bomb_error(const char *message)
{
	fprintf( stderr, "\nFatal Error:\n%s\nProgram halted !!\n\n", message );
	temp_files.Clear();
	exit (1);
} // END void bomb_error

void bomb_error(const char *message, const char *message2)
{
	fprintf( stderr, "\nFatal Error:\n%s %s\nProgram halted !!\n\n", message, message2 );
	temp_files.Clear();
	exit (1);
} // END void bomb_error


void bomb_warning(const char *message)
{
	fprintf( stderr, "\nWarning:\n%s\nNot fatal, but may affect results !!\n\n", message );
} // END void bomb_warning


void bomb_warning(const char *message, const char *message2)
{
	fprintf( stderr, "\nWarning:\n%s %s\nNot fatal, but may affect results !!\n\n", message, message2 );
} // END void bomb_warning

void format_seq(char *seq)
{
	int i, j;
	char c1;
	int len = strlen(seq);

	for (i=0,j=0; i<len; i++) {
		c1 = toupper(seq[i]);
		if ( isalpha(c1) ) seq[j++] = c1;
	}
	seq[j] = 0;
} // END void format_seq

void strrev(char *p)
{
  char *q = p;
  while(q && *q) ++q;
  for(--q; p < q; ++p, --q)
    *p = *p ^ *q,
    *q = *p ^ *q,
    *p = *p ^ *q;
}

/* Quick Sort.
 * Adam Drozdek: Data Structures and Algorithms in C++, 2nd Edition.
 */
void PartialQuickSort( IndexCount *data, int first, int last, int partial )
{
	int lower=first+1, upper=last;
	IndexCount pivot;
	IndexCount val;
	if( first >= last ) return;
	val = data[first];
	data[first] = data[ (first+last)/2 ];
	data[ (first+last)/2 ] = val;
	pivot = data[ first ];

	while( lower <= upper ){
		while( lower <= last && data[lower].count < pivot.count ) lower ++;
		while( pivot.count < data[upper].count ) upper --;
		if( lower < upper ){
			val = data[lower];
			data[lower] = data[upper];
			data[upper] = val;
			upper --;
		}
		lower ++;
	}
	val = data[first];
	data[first] = data[upper];
	data[upper] = val;
	if( first < upper-1 ) PartialQuickSort( data, first, upper-1, partial );
	if( upper >= partial ) return;
	if( upper+1 < last ) PartialQuickSort( data, upper+1, last, partial );
}
Sequence::Sequence()
{
	memset( this, 0, sizeof( Sequence ) );
	distance = 2.0;
}
Sequence::Sequence( const Sequence & other )
{
	//printf( "new: %p  %p\n", this, & other );
	memcpy( this, & other, sizeof( Sequence ) );
	distance = 2.0;
	if( other.data ){
		size = bufsize = other.size;
                size_R2 = 0;
		data = new char[size+1];
		//printf( "data: %p  %p\n", data, other.data );
		data[size] = 0;
		memcpy( data, other.data, size );
		//for (i=0; i<size; i++) data[i] = other.data[i];
	}
	if( other.identifier ){
		int len = strlen( other.identifier );
		identifier = new char[len+1];
		memcpy( identifier, other.identifier, len );
		identifier[len] = 0;
	}
}

// back to back merge for PE
// R1 -> XXXXXXABC ------------------- NMLYYYYYY <--R2
// >R1           >R2
// XXXXXXABC     YYYYYYLMN =====> Merge into
// >R12
// NMLYYYYYYXXXXXXABC
Sequence::Sequence( const Sequence & other, const Sequence & other2, size_t mode )
{
	if (mode != 1) bomb_error("unknown mode");

	//printf( "new: %p  %p\n", this, & other );
	memcpy( this, & other, sizeof( Sequence ) );
	distance = 2.0;

	if( other.data && other2.data ){
		size = bufsize = (other.size + other2.size);
                size_R2 = other2.size;
		data = new char[size+1];
		//printf( "data: %p  %p\n", data, other.data );
		data[size] = 0;     
                data[size_R2] = 0;  
                memcpy( data, other2.data, size_R2); // copy R2 first
                strrev( data );                      // reverse R2 on data
		memcpy( data+size_R2, other.data, size-size_R2 ); // copy R1 to end of R2
		//for (i=0; i<size; i++) data[i] = other.data[i];
		des_begin2 = other2.des_begin;
                tot_length2= other2.tot_length;
	}
        else if ( other.data || other2.data ) {
                bomb_error("Not both PE sequences have data");
        }

	if( other.identifier ){ // only use R1
		int len = strlen( other.identifier );
		identifier = new char[len+1];
		memcpy( identifier, other.identifier, len );
		identifier[len] = 0;
	}
}


Sequence::~Sequence()
{
	//printf( "delete: %p\n", this );
	if( data ) delete[] data;
	if( identifier ) delete[] identifier;
}

void Sequence::Clear()
{
	if( data ) delete[] data;
	/* do not set size to zero here, it is need for writing output */
	bufsize = 0;
	data = NULL;
}

void Sequence::operator=( const char *s )
{
	size = 0; // avoid copying;
	Resize( strlen( s ) );
	strcpy( data, s );
}
void Sequence::operator+=( const char *s )
{
	int m = size; int n = strlen( s );
	Reserve( m + n );
	memcpy( data+m, s, n );
}
void Sequence::Resize( size_t n )
{
	int m = size < n ? size : n;
	size = n;
	if( size != bufsize ){
		char *old = data;
		bufsize = size;
		data = new char[ bufsize + 1 ];
		if ( data == NULL ) bomb_error( "Memory" );
		if ( old ){
			memcpy( data, old, m );
			delete []old;
		}
		if( size ) data[size] = 0;
	}
}
void Sequence::Reserve( size_t n )
{
	int m = size < n ? size : n;
	size = n;
	if( size > bufsize ){
		char *old = data;
		bufsize = size + size/5 + 1;
		data = new char[ bufsize + 1 ];
		if ( data == NULL ) bomb_error( "Memory" );
		if ( old ){
			memcpy( data, old, m );
			delete []old;
		}
	}
	if( size ) data[size] = 0;
}
void Sequence::trim(size_t trim_len) {
    if (trim_len >= size) return;
    size = trim_len;
    if (size) data[size]=0;
}
void Sequence::ConvertBases()
{
	for(size_t i=0; i<size; i++) data[i] = aa2idx[data[i] - 'A'];
}

void Sequence::Swap( Sequence & other )
{
	std::swap(*this, other);
}
int Sequence::Format()
{
	size_t j = 0;
	size_t m = 0;
	while( size && isspace( data[size-1] ) ) size --;
	if( size && data[size-1] == '*' ) size --;
	if( size ) data[size] = 0;
	for (size_t i=0; i<size; i++){
		char ch = data[i];
		m += ! (isalpha( ch ) | isspace( ch ));
	}
	if( m ) return m;
	for (size_t i=0; i<size; i++){
		char ch = data[i];
		if ( isalpha( ch ) ) data[j++] = toupper( ch );
	}
	data[j] = 0;
	size = j;
	return 0;
}

void Sequence::SwapIn()
{
	if( data ) return;
	if( swap == NULL ) bomb_error( "Can not swap in sequence" );
	Resize( size );
	fseek( swap, offset, SEEK_SET );
	if( fread( data, 1, size, swap ) ==0 ) bomb_error( "Can not swap in sequence" );
	data[size] = 0;
}
void Sequence::SwapOut()
{
	if( swap && data ){
		delete[] data;
		bufsize = 0;
		data = NULL;
	}
}
void Sequence::PrintInfo( size_t id, FILE *fout, char *buf )
{
	const auto& tag = options.isEST ? "nt" : "aa";
	bool print = options.print != 0;
	bool strand = options.isEST;
	std::cout << id << "\t" << size << tag << ", >" << (identifier + 1) << "...";
	if( identity ){
		fprintf( fout, " at " );
		if (print) fprintf( fout, "%zu:%zu:%zu:%zu/", coverage[0], coverage[1], coverage[2], coverage[3] );
		if (strand) fprintf( fout, "%c/", (state & IS_MINUS_STRAND) ? '-' : '+' );
		fprintf( fout, "%.2f%%", identity*100 );
		if( options.useDistance ) fprintf( fout, "/%.2f%%", distance*100 );
		fprintf( fout, "\n" );
	}else{
		fprintf( fout, " *\n" );
	}
}

// by liwz gzip version 2019-02
// by liwz
// disable swap option
// change des_begin, des_length, des_length2, dat_length => des_begin, tot_length
// where des_begin is the FILE pointer of sequence record start
//       tot_length is the total bytes of sequence record 
void SequenceDB::Readgz( string file )
{
#ifndef NO_ZLIB
    Sequence one;
    Sequence des;
    auto fin = gzopen(file.c_str(), "r");
    char *buffer = NULL;
    char *res = NULL;
    auto option_l = options.min_length;
    if( fin == NULL ) bomb_error( "Failed to open the database file" );
    Clear();
    buffer = new char[ MAX_LINE_SIZE+1 ];

    while (not gzeof( fin ) || one.size) { /* do not break when the last sequence is not handled */
        buffer[0] = '>';
        if ( (res=gzgets(fin, buffer, MAX_LINE_SIZE)) == NULL && one.size == 0) break;
        if( buffer[0] == '+' ){
            int len = strlen( buffer );
            int len2 = len;
            while( len2 && buffer[len2-1] != '\n' ){
                if ( (res=gzgets(fin, buffer, MAX_LINE_SIZE)) == NULL ) break;
                len2 = strlen( buffer );
                len += len2;
            }
            one.tot_length += len;

            // read next line quality score
            if ( (res=gzgets(fin, buffer, MAX_LINE_SIZE)) == NULL ) bomb_error("can not read quality score after");
            len = strlen( buffer );
            len2 = len;
            while( len2 && buffer[len2-1] != '\n' ){
                if ( (res=gzgets(fin, buffer, MAX_LINE_SIZE)) == NULL ) break;
                len2 = strlen( buffer );
                len += len2;
            }
            one.tot_length += len;
        }else if (buffer[0] == '>' || buffer[0] == '@' || (res==NULL && one.size)) {
            if ( one.size ) { // write previous record
                if( one.identifier == NULL || one.Format() ){
                    printf( "Warning: from file \"%s\",\n", file.c_str() );
                    printf( "Discarding invalid sequence or sequence without identifier and description!\n\n" );
                    if( one.identifier ) printf( "%s\n", one.identifier );
                    printf( "%s\n", one.data );
                    one.size = 0;
                }
                one.index = sequences.size();
                if( one.size > option_l ) {
                    if (options.trim_len    > 0) one.trim(options.trim_len);
                    sequences.push_back( new Sequence( one ) );
                }
            }
            one.size = 0;
            one.tot_length = 0;

            int len = strlen( buffer );
            int len2 = len;
            des.size = 0;
            des += buffer;
            while( len2 && buffer[len2-1] != '\n' ){
                if ( (res=gzgets(fin, buffer, MAX_LINE_SIZE)) == NULL ) break;
                des += buffer;
                len2 = strlen( buffer );
                len += len2;
            }
            size_t offset = gztell( fin );
            one.des_begin = offset - len;
            one.tot_length += len;              // count first line

            size_t i = 0;
            if( des.data[i] == '>' || des.data[i] == '@' || des.data[i] == '+' ) i += 1;
            if( des.data[i] == ' ' or des.data[i] == '\t' ) i += 1;
            if( options.des_len and options.des_len < des.size ) des.size = options.des_len;
            while( i < des.size and ! isspace( des.data[i] ) ) i += 1;
            des.data[i] = 0;
            one.identifier = des.data;
        } else {
            one.tot_length += strlen(buffer);  one += buffer;
        }
    }
    one.identifier = NULL;
    delete[] buffer;
    gzclose( fin );
#else
    bomb_error("this program was not compiled with zlib");
#endif
}



// by liwz
// disable swap option
// change des_begin, des_length, des_length2, dat_length => des_begin, tot_length
// where des_begin is the FILE pointer of sequence record start
//       tot_length is the total bytes of sequence record 
void SequenceDB::Read( const char *sfile )
{
    const string file{sfile};
    if (file.ends_with(".gz")) {
        Readgz(file);
        return;
    }

    Sequence one;
    Sequence des;
    FILE *fin = fopen( file.c_str(), "rb" );
    char *buffer = NULL;
    char *res = NULL;
    if( fin == NULL ) bomb_error( "Failed to open the database file" );
    Clear();
    buffer = new char[ MAX_LINE_SIZE+1 ];

    while (not feof( fin ) || one.size) { /* do not break when the last sequence is not handled */
        buffer[0] = '>';
        if ( (res=fgets( buffer, MAX_LINE_SIZE, fin )) == NULL && one.size == 0) break;
        if( buffer[0] == '+' ){
            int len = strlen( buffer );
            int len2 = len;
            while( len2 && buffer[len2-1] != '\n' ){
                if ( (res=fgets( buffer, MAX_LINE_SIZE, fin )) == NULL ) break;
                len2 = strlen( buffer );
                len += len2;
            }
            one.tot_length += len;

            // read next line quality score
            if ( (res=fgets( buffer, MAX_LINE_SIZE, fin )) == NULL ) bomb_error("can not read quality score after");
            len = strlen( buffer );
            len2 = len;
            while( len2 && buffer[len2-1] != '\n' ){
                if ( (res=fgets( buffer, MAX_LINE_SIZE, fin )) == NULL ) break;
                len2 = strlen( buffer );
                len += len2;
            }
            one.tot_length += len;
        }else if (buffer[0] == '>' || buffer[0] == '@' || (res==NULL && one.size)) {
            if ( one.size ) { // write previous record
                if( one.identifier == NULL || one.Format() ){
                    cout << "Warning: from file \"" << file << "\"," << endl;
                    printf( "Discarding invalid sequence or sequence without identifier and description!\n\n" );
                    if( one.identifier ) cout << one.identifier << endl;
                    cout << one.data << endl;
                    one.size = 0;
                }
                one.index = sequences.size();
                if( one.size > options.min_length ) {
                    if (options.trim_len    > 0) one.trim(options.trim_len);
                    sequences.Append( new Sequence( one ) ); 
                }
            }
            one.size = 0;
            one.tot_length = 0;

            int len = strlen( buffer );
            int len2 = len;
            des.size = 0;
            des += buffer;
            while( len2 && buffer[len2-1] != '\n' ){
                if ( (res=fgets( buffer, MAX_LINE_SIZE, fin )) == NULL ) break;
                des += buffer;
                len2 = strlen( buffer );
                len += len2;
            }
            size_t offset = ftell( fin );
            one.des_begin = offset - len;
            one.tot_length += len;              // count first line

            size_t i = 0;
            if( des.data[i] == '>' || des.data[i] == '@' || des.data[i] == '+' ) i += 1;
            if( des.data[i] == ' ' or des.data[i] == '\t' ) i += 1;
            if( options.des_len and options.des_len < des.size ) des.size = options.des_len;
            while( i < des.size and ! isspace( des.data[i] ) ) i += 1;
            des.data[i] = 0;
            one.identifier = des.data;
        } else {
            one.tot_length += strlen(buffer);  one += buffer;
        }
    }
#if 0
    int i, n = 0;
    for(i=0; i<sequences.size(); i++) n += sequences[i].bufsize + 4;
    cout<<n<<"\t"<<sequences.capacity() * sizeof(Sequence)<<endl;
    int i;
    scanf( "%i", & i );
#endif
    one.identifier = NULL;
    delete[] buffer;
    fclose( fin );
}


// by liwz gzip version 2019-02
// PE reads liwz, disable swap option
void SequenceDB::Readgz(string file, string file2 )
{
#ifndef NO_ZLIB
    Sequence one, two;
    Sequence des;
    gzFile fin = gzopen(file.c_str(), "r");
    gzFile fin2= gzopen(file2.c_str(),"r");
    char *buffer = NULL;
    char *buffer2= NULL;
    char *res = NULL;
    char *res2= NULL;
    size_t option_l = options.min_length;
    if( fin == NULL ) bomb_error( "Failed to open the database file" );
    if( fin2== NULL ) bomb_error( "Failed to open the database file" );
    Clear();
    buffer = new char[ MAX_LINE_SIZE+1 ];
    buffer2= new char[ MAX_LINE_SIZE+1 ];

    while (((not gzeof( fin )) && (not gzeof( fin2)) ) || (one.size && two.size)) { /* do not break when the last sequence is not handled */
        buffer[0] = '>'; res =gzgets(fin,  buffer,  MAX_LINE_SIZE);
        buffer2[0]= '>'; res2=gzgets(fin2, buffer2, MAX_LINE_SIZE);

        if ( (res      == NULL) && (res2     != NULL)) bomb_error( "Paired input files have different number sequences" );
        if ( (res      != NULL) && (res2     == NULL)) bomb_error( "Paired input files have different number sequences" );
        if ( (one.size == 0   ) && (two.size >     0)) bomb_error( "Paired input files have different number sequences" );
        if ( (one.size >  0   ) && (two.size ==    0)) bomb_error( "Paired input files have different number sequences" );
        if ( (res      == NULL) && (one.size ==    0)) break;

        if( buffer[0] == '+' ){ // fastq 3rd line
            // file 1
            int len = strlen( buffer ); 
            int len2 = len;
            while( len2 && buffer[len2-1] != '\n' ){ // read until the end of the line
                if ( (res=gzgets(fin, buffer, MAX_LINE_SIZE)) == NULL ) break;
                len2 = strlen( buffer );
                len += len2;
            }
            one.tot_length += len;

            // read next line quality score
            if ( (res=gzgets(fin, buffer, MAX_LINE_SIZE)) == NULL ) bomb_error("can not read quality score after");
            len = strlen( buffer );
            len2 = len;
            while( len2 && buffer[len2-1] != '\n' ){
                if ( (res=gzgets(fin, buffer, MAX_LINE_SIZE)) == NULL ) break;
                len2 = strlen( buffer );
                len += len2;
            }
            one.tot_length += len;

            // file 2
            len = strlen( buffer2 );
            len2 = len;
            while( len2 && buffer2[len2-1] != '\n' ){ // read until the end of the line
                if ( (res2=gzgets(fin2, buffer2, MAX_LINE_SIZE)) == NULL ) break;
                len2 = strlen( buffer2 );
                len += len2;
            }
            two.tot_length += len;

            // read next line quality score
            if ( (res2=gzgets(fin2, buffer2, MAX_LINE_SIZE)) == NULL ) bomb_error("can not read quality score after");
            len = strlen( buffer2 );
            len2 = len;
            while( len2 && buffer2[len2-1] != '\n' ){
                if ( (res2=gzgets(fin2, buffer2, MAX_LINE_SIZE)) == NULL ) break;
                len2 = strlen( buffer2 );
                len += len2;
            }
            two.tot_length += len;

        }else if (buffer[0] == '>' || buffer[0] == '@' || (res==NULL && one.size)) {
            if ( one.size && two.size ) { // write previous record
                if( one.identifier == NULL || one.Format() ){
                    printf( "Warning: from file \"%s\",\n", file.c_str() );
                    printf( "Discarding invalid sequence or sequence without identifier and description!\n\n" );
                    if( one.identifier ) printf( "%s\n", one.identifier );
                    printf( "%s\n", one.data );
                    one.size=0; two.size=0;
                }
                if( two.identifier == NULL || two.Format() ){
                    printf( "Warning: from file \"%s\",\n", file2.c_str() );
                    printf( "Discarding invalid sequence or sequence without identifier and description!\n\n" );
                    if( two.identifier ) printf( "%s\n", two.identifier );
                    printf( "%s\n", two.data );
                    one.size=0; two.size = 0;
                }
                one.index = sequences.size();
                if( (one.size + two.size)> option_l ) {
                    if (options.trim_len    > 0) one.trim(options.trim_len);
                    if (options.trim_len_R2 > 0) two.trim(options.trim_len_R2);
                    sequences.Append( new Sequence( one, two, 1 ) ); 
                }
            }
            // R1
            one.size = 0;
            one.tot_length = 0;

            int len = strlen( buffer );
            int len2 = len;
            des.size = 0;
            des += buffer;
            while( len2 && buffer[len2-1] != '\n' ){
                if ( (res=gzgets(fin, buffer, MAX_LINE_SIZE)) == NULL ) break;
                des += buffer;
                len2 = strlen( buffer );
                len += len2;
            }
            size_t offset = gztell( fin );    
            one.des_begin = offset - len; // offset of ">" or "@" 
            one.tot_length += len;              // count first line

            size_t i = 0;
            if( des.data[i] == '>' || des.data[i] == '@' || des.data[i] == '+' ) i += 1;
            if( des.data[i] == ' ' or des.data[i] == '\t' ) i += 1;
            if( options.des_len and options.des_len < des.size ) des.size = options.des_len;
            while( i < des.size and ! isspace( des.data[i] ) ) i += 1;
            des.data[i] = 0;                   // find first non-space letter
            one.identifier = des.data;

            // R2
            two.size = 0;
            two.tot_length = 0;

            len = strlen( buffer2 );
            len2 = len;
            while( len2 && buffer2[len2-1] != '\n' ){
                if ( (res=gzgets(fin2, buffer2, MAX_LINE_SIZE)) == NULL ) break;
                len2 = strlen( buffer2 );
                len += len2;
            }
            offset = gztell( fin2 );
            two.des_begin = offset - len;
            two.tot_length += len;              // count first line
            two.identifier = des.data;
        } else {
            one.tot_length += strlen(buffer);  one += buffer;
            two.tot_length+= strlen(buffer2); two+= buffer2;
        }
    }
    one.identifier = NULL;
    two.identifier = NULL;
    delete[] buffer;
    gzclose( fin );
    delete[] buffer2;
    gzclose( fin2 );
#else
    bomb_error("this program was not compiled with zlib");
#endif

}

// PE reads liwz, disable swap option
void SequenceDB::Read( const char *file, const char *file2 )
{
    int f_len = strlen(file);
    int f_len2= strlen(file2);
    if (strcmp(file + f_len - 3, ".gz") == 0 ) {
        if ( strcmp(file2 + f_len2 - 3, ".gz") ) bomb_error( "Both input files need to be in .gz format" );
        Readgz(file, file2);
        return;
    }

    Sequence one, two;
    Sequence des;
    FILE *fin = fopen( file, "rb" );
    FILE *fin2= fopen( file2,"rb" );
    char *buffer = NULL;
    char *buffer2= NULL;
    char *res = NULL;
    char *res2= NULL;
    size_t option_l = options.min_length;
    if( fin == NULL ) bomb_error( "Failed to open the database file" );
    if( fin2== NULL ) bomb_error( "Failed to open the database file" );
    Clear();
    buffer = new char[ MAX_LINE_SIZE+1 ];
    buffer2= new char[ MAX_LINE_SIZE+1 ];

    while (((not feof( fin )) && (not feof( fin2)) ) || (one.size && two.size)) { /* do not break when the last sequence is not handled */
        buffer[0] = '>'; res =fgets( buffer,  MAX_LINE_SIZE, fin  );
        buffer2[0]= '>'; res2=fgets( buffer2, MAX_LINE_SIZE, fin2 );

        if ( (res      == NULL) && (res2     != NULL)) bomb_error( "Paired input files have different number sequences" );
        if ( (res      != NULL) && (res2     == NULL)) bomb_error( "Paired input files have different number sequences" );
        if ( (one.size == 0   ) && (two.size >     0)) bomb_error( "Paired input files have different number sequences" );
        if ( (one.size >  0   ) && (two.size ==    0)) bomb_error( "Paired input files have different number sequences" );
        if ( (res      == NULL) && (one.size ==    0)) break;

        if( buffer[0] == '+' ){ // fastq 3rd line
            // file 1
            int len = strlen( buffer ); 
            int len2 = len;
            while( len2 && buffer[len2-1] != '\n' ){ // read until the end of the line
                if ( (res=fgets( buffer, MAX_LINE_SIZE, fin )) == NULL ) break;
                len2 = strlen( buffer );
                len += len2;
            }
            one.tot_length += len;

            // read next line quality score
            if ( (res=fgets( buffer, MAX_LINE_SIZE, fin )) == NULL ) bomb_error("can not read quality score after");
            len = strlen( buffer );
            len2 = len;
            while( len2 && buffer[len2-1] != '\n' ){
                if ( (res=fgets( buffer, MAX_LINE_SIZE, fin )) == NULL ) break;
                len2 = strlen( buffer );
                len += len2;
            }
            one.tot_length += len;

            // file 2
            len = strlen( buffer2 );
            len2 = len;
            while( len2 && buffer2[len2-1] != '\n' ){ // read until the end of the line
                if ( (res2=fgets( buffer2, MAX_LINE_SIZE, fin2 )) == NULL ) break;
                len2 = strlen( buffer2 );
                len += len2;
            }
            two.tot_length += len;

            // read next line quality score
            if ( (res2=fgets( buffer2, MAX_LINE_SIZE, fin2 )) == NULL ) bomb_error("can not read quality score after");
            len = strlen( buffer2 );
            len2 = len;
            while( len2 && buffer2[len2-1] != '\n' ){
                if ( (res2=fgets( buffer2, MAX_LINE_SIZE, fin2 )) == NULL ) break;
                len2 = strlen( buffer2 );
                len += len2;
            }
            two.tot_length += len;

        }else if (buffer[0] == '>' || buffer[0] == '@' || (res==NULL && one.size)) {
            if ( one.size && two.size ) { // write previous record
                if( one.identifier == NULL || one.Format() ){
                    printf( "Warning: from file \"%s\",\n", file );
                    printf( "Discarding invalid sequence or sequence without identifier and description!\n\n" );
                    if( one.identifier ) printf( "%s\n", one.identifier );
                    printf( "%s\n", one.data );
                    one.size=0; two.size=0;
                }
                if( two.identifier == NULL || two.Format() ){
                    printf( "Warning: from file \"%s\",\n", file2 );
                    printf( "Discarding invalid sequence or sequence without identifier and description!\n\n" );
                    if( two.identifier ) printf( "%s\n", two.identifier );
                    printf( "%s\n", two.data );
                    one.size=0; two.size = 0;
                }
                one.index = sequences.size();
                if( (one.size + two.size)> option_l ) {
                    if (options.trim_len    > 0) one.trim(options.trim_len);
                    if (options.trim_len_R2 > 0) two.trim(options.trim_len_R2);
                    sequences.Append( new Sequence( one, two, 1 ) ); 
                }
            }
            // R1
            one.size = 0;
            one.tot_length = 0;

            size_t len = strlen( buffer );
            size_t len2 = len;
            des.size = 0;
            des += buffer;
            while( len2 && buffer[len2-1] != '\n' ){
                if ( (res=fgets( buffer, MAX_LINE_SIZE, fin )) == NULL ) break;
                des += buffer;
                len2 = strlen( buffer );
                len += len2;
            }
            size_t offset = ftell( fin );    
            one.des_begin = offset - len; // offset of ">" or "@" 
            one.tot_length += len;              // count first line

            size_t i = 0;
            if( des.data[i] == '>' || des.data[i] == '@' || des.data[i] == '+' ) i += 1;
            if( des.data[i] == ' ' or des.data[i] == '\t' ) i += 1;
            if( options.des_len and options.des_len < des.size ) des.size = options.des_len;
            while( i < des.size and ! isspace( des.data[i] ) ) i += 1;
            des.data[i] = 0;                   // find first non-space letter
            one.identifier = des.data;

            // R2
            two.size = 0;
            two.tot_length = 0;

            len = strlen( buffer2 );
            len2 = len;
            while( len2 && buffer2[len2-1] != '\n' ){
                if ( (res=fgets( buffer2, MAX_LINE_SIZE, fin2 )) == NULL ) break;
                len2 = strlen( buffer2 );
                len += len2;
            }
            offset = ftell( fin2 );
            two.des_begin = offset - len;
            two.tot_length += len;              // count first line
            two.identifier = des.data;
        } else {
            one.tot_length += strlen(buffer);  one += buffer;
            two.tot_length+= strlen(buffer2); two+= buffer2;
        }
    }
#if 0
    int i, n = 0;
    for(i=0; i<sequences.size(); i++) n += sequences[i].bufsize + 4;
    cout<<n<<"\t"<<sequences.capacity() * sizeof(Sequence)<<endl;
    int i;
    scanf( "%i", & i );
#endif
    one.identifier = NULL;
    two.identifier = NULL;
    delete[] buffer;
    fclose( fin );
    delete[] buffer2;
    fclose( fin2 );
}

#if 0
void SequenceDB::Sort( int first, int last )
{
	int lower=first+1, upper=last;
	Sequence *pivot;
	Sequence *val;
	if( first >= last ) return;
	val = sequences[first];
	sequences[first] = sequences[ (first+last)/2 ];
	sequences[ (first+last)/2 ] = val;
	pivot = sequences[ first ];

	while( lower <= upper ){
		while( lower <= last && sequences[lower]->stats < pivot->stats ) lower ++;
		while( pivot->stats < sequences[upper]->stats ) upper --;
		if( lower < upper ){
			val = sequences[lower];
			sequences[lower] = sequences[upper];
			sequences[upper] = val;
			upper --;
		}
		lower ++;
	}
	val = sequences[first];
	sequences[first] = sequences[upper];
	sequences[upper] = val;
	if( first < upper-1 ) Sort( first, upper-1 );
	if( upper+1 < last ) Sort( upper+1, last );
}
#endif
void SequenceDB::SortDivide( Options & options, bool sort )
{
	int i;
	unsigned int len;
	int N = sequences.size();
	total_letter=0;
	total_desc=0;
	max_len = 0;
	min_len = (size_t)-1;
	for (i=0; i<N; i++) {
		Sequence *seq = sequences[i];
		len = seq->size;
		total_letter += len;
		if (len > max_len) max_len = len;
		if (len < min_len) min_len = len;
		if (seq->swap == NULL) seq->ConvertBases();
		if( seq->identifier ) total_desc += strlen( seq->identifier );
	}
	options.max_entries = max_len * MAX_TABLE_SEQ;
	if (max_len >= 65536 and sizeof(INTs) <=2) 
		bomb_warning("Some seqs longer than 65536, you may define LONG_SEQ");

	if (max_len > MAX_SEQ ) 
		bomb_warning("Some seqs are too long, please rebuild the program with make parameter "
				"MAX_SEQ=new-maximum-length (e.g. make MAX_SEQ=10000000)");

	cout << "longest and shortest : " << max_len << " and " << min_len << endl;
	cout << "Total letters: " << total_letter << endl;
	// END change all the NR_seq to iseq

	len_n50 = (max_len + min_len) / 2; // will be properly set, if sort is true;
	if( sort ){
		// **************************** Form NR_idx[], Sort them from Long to short
		long long sum = 0;
		int M = max_len - min_len + 1;
		std::vector<int> count( M, 0 ); // count for each size = max_len - i
		std::vector<int> accum( M, 0 ); // count for all size > max_len - i
		std::vector<int> offset( M, 0 ); // offset from accum[i] when filling sorting
		std::vector<Sequence*> sorting( N ); // TODO: use a smaller class if this consumes to much memory!

		for (i=0; i<N; i++) count[ max_len - sequences[i]->size ] ++;
		for (i=1; i<M; i++) accum[i] = accum[i-1] + count[i-1];
		for (i=0; i<M; i++){
			sum += (max_len - i) * count[i];
			if( sum >= (total_letter>>1) ){
				len_n50 = max_len - i;
				break;
			}
		}
		for (i=0; i<N; i++){
			int len = max_len - sequences[i]->size;
			int id = accum[len] + offset[len];
			//sequences[i].index = id;
			sorting[id] = sequences[i];
			offset[len] ++;
		}
		options.max_entries = 0;
		for (i=0; i<N; i++){
			sequences[i] = sorting[i];
			if( i < MAX_TABLE_SEQ ) options.max_entries += sequences[i]->size;
		}
#if 0
		if( options.isEST ){
			int start = 0;
			for (i=0; i<M; i++){
				Sort( start, accum[i] );
				start = accum[i];
			}
		}
#endif
		cout << "Sequences have been sorted" << endl;
		// END sort them from long to short
	}
}// END sort_seqs_divide_segs

void SequenceDB::DivideSave( const char *db, const char *newdb, int n )
{
	if( n == 0 or sequences.size() ==0 ) return;

	size_t max_seg = total_letter / n + sequences[0]->size;
	if( max_seg >= MAX_BIN_SWAP ) max_seg = (size_t) MAX_BIN_SWAP;

	FILE *fin = fopen( db, "rb" );
	char *buf = new char[MAX_LINE_SIZE+1];
	char outfile[512];
	size_t seg_size = 0;
	int i, j, count, rest, seg = 0;
	sprintf( outfile, "%s-%i", newdb, 0 );
	FILE *fout = fopen( outfile, "w+" );
	n = sequences.size();
	for (i=0; i<n; i++){
		Sequence *seq = sequences[i];
		fseek( fin, seq->des_begin, SEEK_SET );

		seg_size += seq->size;
		if( seg_size >= max_seg ){
			seg += 1;
			sprintf( outfile, "%s-%i", newdb, seg );
			fclose( fout );
			fout = fopen( outfile, "w+" );
			seg_size = seq->size;
		}

		count = seq->tot_length / MAX_LINE_SIZE;
		rest  = seq->tot_length % MAX_LINE_SIZE;
		//printf( "count = %6i,  rest = %6i\n", count, rest );
		for (j=0; j<count; j++){
			if( fread( buf, 1, MAX_LINE_SIZE, fin ) ==0 ) bomb_error( "Can not swap in sequence" );
			fwrite( buf, 1, MAX_LINE_SIZE, fout );
		}
		if( rest ){
			if( fread( buf, 1, rest, fin ) ==0 ) bomb_error( "Can not swap in sequence" );
			fwrite( buf, 1, rest, fout );
		}
	}
	fclose( fin );
	fclose( fout );
	delete []buf;
}

// input db is gzipped
void SequenceDB::WriteClustersgz( const char *db, const char *newdb )
{
#ifndef NO_ZLIB
    gzFile fin = gzopen(db, "r");
	FILE *fout = fopen( newdb, "w+" );
	int i, j, n = rep_seqs.size();
	int count, rest;
	char *buf = new char[MAX_LINE_SIZE+1];
	vector<uint64_t> sorting( n );
	if( fin == NULL || fout == NULL ) bomb_error( "file opening failed" );
	for (i=0; i<n; i++) sorting[i] = ((uint64_t)sequences[ rep_seqs[i] ]->index << 32) | rep_seqs[i];
	std::sort( sorting.begin(), sorting.end() );
	for (i=0; i<n; i++){
		Sequence *seq = sequences[ sorting[i] & 0xffffffff ];
		gzseek( fin, seq->des_begin, SEEK_SET );

		count = seq->tot_length / MAX_LINE_SIZE;
		rest  = seq->tot_length % MAX_LINE_SIZE;
		//printf( "count = %6i,  rest = %6i\n", count, rest );
		for (j=0; j<count; j++){
			if( gzread(fin, buf, MAX_LINE_SIZE) ==0 ) bomb_error( "Can not swap in sequence" );
			fwrite( buf, 1, MAX_LINE_SIZE, fout );
		}
		if( rest ){
			if( gzread(fin, buf, rest) ==0 ) bomb_error( "Can not swap in sequence" );
			fwrite( buf, 1, rest, fout );
		}
	}
	gzclose( fin );
	fclose( fout );
	delete []buf;
#else
    bomb_error("this program was not compiled with zlib");
#endif

}

void SequenceDB::WriteClusters( const char *db, const char *newdb )
{
    int f_len = strlen(db);
    if (strcmp(db + f_len - 3, ".gz") == 0 ) {
        WriteClustersgz(db, newdb);
        return;
    }

	FILE *fin = fopen( db, "rb" );
	FILE *fout = fopen( newdb, "w+" );
	int i, j, n = rep_seqs.size();
	int count, rest;
	char *buf = new char[MAX_LINE_SIZE+1];
	vector<uint64_t> sorting( n );
	if( fin == NULL || fout == NULL ) bomb_error( "file opening failed" );
	for (i=0; i<n; i++) sorting[i] = ((uint64_t)sequences[ rep_seqs[i] ]->index << 32) | rep_seqs[i];
	std::sort( sorting.begin(), sorting.end() );
	for (i=0; i<n; i++){
		Sequence *seq = sequences[ sorting[i] & 0xffffffff ];
		fseek( fin, seq->des_begin, SEEK_SET );

		count = seq->tot_length / MAX_LINE_SIZE;
		rest  = seq->tot_length % MAX_LINE_SIZE;
		//printf( "count = %6i,  rest = %6i\n", count, rest );
		for (j=0; j<count; j++){
			if( fread( buf, 1, MAX_LINE_SIZE, fin ) ==0 ) bomb_error( "Can not swap in sequence" );
			fwrite( buf, 1, MAX_LINE_SIZE, fout );
		}
		if( rest ){
			if( fread( buf, 1, rest, fin ) ==0 ) bomb_error( "Can not swap in sequence" );
			fwrite( buf, 1, rest, fout );
		}
	}
	fclose( fin );
	fclose( fout );
	delete []buf;
}


// input db is gzipped
// liwz PE output
void SequenceDB::WriteClustersgz( const char *db, const char *db_pe, const char *newdb, const char *newdb_pe )
{
#ifndef NO_ZLIB
    gzFile fin    = gzopen(db,    "r");
	gzFile fin_pe = gzopen(db_pe, "r");
	FILE *fout = fopen( newdb, "w+" );
	FILE *fout_pe = fopen( newdb_pe, "w+" );
	int i, j, n = rep_seqs.size();
	int count, rest;
	char *buf = new char[MAX_LINE_SIZE+1];
	vector<uint64_t> sorting( n );
	if( fin == NULL || fout == NULL ) bomb_error( "file opening failed" );
	if( fin_pe == NULL || fout_pe == NULL ) bomb_error( "file opening failed" );
	for (i=0; i<n; i++) sorting[i] = ((uint64_t)sequences[ rep_seqs[i] ]->index << 32) | rep_seqs[i];
	std::sort( sorting.begin(), sorting.end() );

        //sort fasta / fastq
        int *clstr_size;
        int *clstr_idx1;
        if (options.sort_outputf) {
            clstr_size = new int[n];
            clstr_idx1 = new int[n];
            for (i=0; i<n; i++) { 
                clstr_size[i] = 0;
                clstr_idx1[i]  = i;
            }

            int N = sequences.size();
            for (i=0; i<N; i++) { 
                int id = sequences[i]->cluster_id;
                if (id < 0) continue;
                if (id >=n) continue;
                clstr_size[id]++;
            }
            quick_sort_idxr(clstr_size, clstr_idx1, 0, n-1);
        }

	for (i=0; i<n; i++){
		Sequence *seq = sequences[ sorting[i] & 0xffffffff ];
                if (options.sort_outputf) seq = sequences[  rep_seqs[ clstr_idx1[i] ] ];
                //R1
		gzseek( fin, seq->des_begin, SEEK_SET );

		count = seq->tot_length / MAX_LINE_SIZE;
		rest  = seq->tot_length % MAX_LINE_SIZE;
		//printf( "count = %6i,  rest = %6i\n", count, rest );
		for (j=0; j<count; j++){
			if( gzread(fin, buf, MAX_LINE_SIZE) ==0 ) bomb_error( "Can not swap in sequence" );
			fwrite( buf, 1, MAX_LINE_SIZE, fout );
		}
		if( rest ){
			if( gzread(fin, buf, rest) ==0 ) bomb_error( "Can not swap in sequence" );
			fwrite( buf, 1, rest, fout );
		}

                //R2
		gzseek( fin_pe, seq->des_begin2, SEEK_SET );

		count = seq->tot_length2 / MAX_LINE_SIZE;
		rest  = seq->tot_length2 % MAX_LINE_SIZE;
		//printf( "count = %6i,  rest = %6i\n", count, rest );
		for (j=0; j<count; j++){
			if( gzread(fin_pe, buf, MAX_LINE_SIZE) ==0 ) bomb_error( "Can not swap in sequence" );
			fwrite( buf, 1, MAX_LINE_SIZE, fout_pe );
		}
		if( rest ){
			if( gzread(fin_pe, buf, rest) ==0 ) bomb_error( "Can not swap in sequence" );
			fwrite( buf, 1, rest, fout_pe );
		}

	}
	gzclose( fin );
	gzclose( fin_pe );
	fclose( fout );
	fclose( fout_pe );
	delete []buf;
#else
    bomb_error("this program was not compiled with zlib");
#endif
}


// liwz PE output
void SequenceDB::WriteClusters( const char *db, const char *db_pe, const char *newdb, const char *newdb_pe )
{
    int f_len = strlen(db); 
    if (strcmp(db + f_len - 3, ".gz") == 0 ) {
        WriteClustersgz(db, db_pe, newdb, newdb_pe);
        return;
    }

	FILE *fin = fopen( db, "rb" );
	FILE *fout = fopen( newdb, "w+" );
	FILE *fin_pe = fopen( db_pe, "rb" );
	FILE *fout_pe = fopen( newdb_pe, "w+" );
	int i, j, n = rep_seqs.size();
	int count, rest;
	char *buf = new char[MAX_LINE_SIZE+1];
	vector<uint64_t> sorting( n );
	if( fin == NULL || fout == NULL ) bomb_error( "file opening failed" );
	if( fin_pe == NULL || fout_pe == NULL ) bomb_error( "file opening failed" );
	for (i=0; i<n; i++) sorting[i] = ((uint64_t)sequences[ rep_seqs[i] ]->index << 32) | rep_seqs[i];
	std::sort( sorting.begin(), sorting.end() );

        //sort fasta / fastq
        int *clstr_size;
        int *clstr_idx1;
        if (options.sort_outputf) {
            clstr_size = new int[n];
            clstr_idx1 = new int[n];
            for (i=0; i<n; i++) { 
                clstr_size[i] = 0;
                clstr_idx1[i]  = i;
            }

            int N = sequences.size();
            for (i=0; i<N; i++) { 
                int id = sequences[i]->cluster_id;
                if (id < 0) continue;
                if (id >=n) continue;
                clstr_size[id]++;
            }
            quick_sort_idxr(clstr_size, clstr_idx1, 0, n-1);
        }

	for (i=0; i<n; i++){
		Sequence *seq = sequences[ sorting[i] & 0xffffffff ];
                if (options.sort_outputf) seq = sequences[  rep_seqs[ clstr_idx1[i] ] ];
                //R1
		fseek( fin, seq->des_begin, SEEK_SET );

		count = seq->tot_length / MAX_LINE_SIZE;
		rest  = seq->tot_length % MAX_LINE_SIZE;
		//printf( "count = %6i,  rest = %6i\n", count, rest );
		for (j=0; j<count; j++){
			if( fread( buf, 1, MAX_LINE_SIZE, fin ) ==0 ) bomb_error( "Can not swap in sequence" );
			fwrite( buf, 1, MAX_LINE_SIZE, fout );
		}
		if( rest ){
			if( fread( buf, 1, rest, fin ) ==0 ) bomb_error( "Can not swap in sequence" );
			fwrite( buf, 1, rest, fout );
		}

                //R2
		fseek( fin_pe, seq->des_begin2, SEEK_SET );

		count = seq->tot_length2 / MAX_LINE_SIZE;
		rest  = seq->tot_length2 % MAX_LINE_SIZE;
		//printf( "count = %6i,  rest = %6i\n", count, rest );
		for (j=0; j<count; j++){
			if( fread( buf, 1, MAX_LINE_SIZE, fin_pe ) ==0 ) bomb_error( "Can not swap in sequence" );
			fwrite( buf, 1, MAX_LINE_SIZE, fout_pe );
		}
		if( rest ){
			if( fread( buf, 1, rest, fin_pe ) ==0 ) bomb_error( "Can not swap in sequence" );
			fwrite( buf, 1, rest, fout_pe );
		}

	}
	fclose( fin );
	fclose( fout );
	fclose( fin_pe );
	fclose( fout_pe );
	delete []buf;
}

void SequenceDB::WriteExtra1D( const Options & options )
{
	string db_clstr = options.output + ".clstr";
	string db_clstr_bak = options.output + ".bak.clstr";
	int i, i0, k, N = sequences.size();
	vector<long long> sorting( N );
	for (i=0; i<N; i++) sorting[i] = ((long long)sequences[i]->index << 32) | i;
	std::sort( sorting.begin(), sorting.end() );

	FILE *fout;
	char *buf = new char[ MAX_DES + 1 ];

	if( options.backupFile ){
		fout = fopen( db_clstr_bak.c_str(), "w+" );
		for (i=0; i<N; i++) {
			Sequence *seq = sequences[ sorting[i] & 0xffffffff ];
			seq->PrintInfo( seq->cluster_id, fout, buf );
		}
		fclose( fout );
	}

	cout << "writing clustering information" << endl;
	int M = rep_seqs.size();
	std::vector<std::vector<int> > clusters( M );
	for (i=0; i<N; i++){
		int k = sorting[i] & 0xffffffff;
		int id = sequences[k]->cluster_id;
		clusters[id].Append( k );
	}

	fout = fopen( db_clstr.c_str(), "w+" );

        if (options.sort_output) {
            int *clstr_size = new int[M];
            int *clstr_idx1 = new int[M];

            for (i=0; i<M; i++) { 
                clstr_size[i] = (int)clusters[i].size();
                clstr_idx1[i]  = i;
            }
            quick_sort_idxr(clstr_size, clstr_idx1, 0, M-1);

  	    for (i=0; i<M; i++) {
                i0 = clstr_idx1[i];
		fprintf( fout, ">Cluster %i\n", i );
		for (k=0; k<(int)clusters[i0].size(); k++)
			sequences[ clusters[i0][k] ]->PrintInfo( k, fout, buf );
	    }   
        }
        else {
  	    for (i=0; i<M; i++) {
		fprintf( fout, ">Cluster %i\n", i );
		for (k=0; k<(int)clusters[i].size(); k++)
			sequences[ clusters[i][k] ]->PrintInfo( k, fout, buf );
	    }   

        }

	delete []buf;
}
void SequenceDB::WriteExtra2D( SequenceDB & other)
{
	string db_clstr = options.output + ".clstr";
	string db_clstr_bak = options.output + ".bak.clstr";
	int i, k, N = other.sequences.size();
	int N2 = sequences.size();
	vector<long long> sorting( N );
	for (i=0; i<N; i++) sorting[i] = ((long long)other.sequences[i]->index << 32) | i;
	std::sort( sorting.begin(), sorting.end() );

	FILE *fout;
	char *buf = new char[ MAX_DES + 1 ];
	if( options.backupFile ){
		fout = fopen( db_clstr_bak.c_str(), "w+" );
		for (i=0; i<N; i++) {
			Sequence *seq = other.sequences[ sorting[i] & 0xffffffff ];
			seq->PrintInfo( seq->cluster_id, fout, buf );
		}
		for (i=0; i<N2; i++) {
			Sequence *seq = sequences[i];
			if( seq->state & IS_REDUNDANT ) seq->PrintInfo( seq->cluster_id, fout, buf );
		}
		fclose( fout );
	}

	cout << "writing clustering information" << endl;
	std::vector<std::vector<int> > clusters( N );
	for (i=0; i<N2; i++){
		int id = sequences[i]->cluster_id;
		if( sequences[i]->state & IS_REDUNDANT ) clusters[id].Append( i );
	}

	fout = fopen( db_clstr.c_str(), "w+" );
	for (i=0; i<N; i++) {
		Sequence *seq = other.sequences[ i ];
		fprintf( fout, ">Cluster %i\n", i );
		seq->PrintInfo( 0, fout, buf );
		for (k=0; k<(int)clusters[i].size(); k++)
			sequences[ clusters[i][k] ]->PrintInfo( k+1, fout, buf );
	}
	delete []buf;
}


void WorkingParam::ControlShortCoverage( int len)
{
	len_eff = len;
	aln_cover_flag = 0;
	if ((options.short_coverage > 0.0) || (options.min_control>0) ) { // has alignment coverage control
		aln_cover_flag = 1;
		min_aln_lenS = (int) (double(len) * options.short_coverage);
		if ( len-options.short_control > min_aln_lenS) min_aln_lenS = len-options.short_control;
		if ( options.min_control > min_aln_lenS) min_aln_lenS = options.min_control;
	}
	if (options.global_identity == 0) len_eff = min_aln_lenS; //global_identity==0
}
void WorkingParam::ControlLongCoverage( int len2)
{
	if (aln_cover_flag) {
		min_aln_lenL = (int) (double(len2) * options.long_coverage);
		if ( len2-options.long_control > min_aln_lenL) min_aln_lenL = len2-options.long_control;
		if ( options.min_control > min_aln_lenL) min_aln_lenL = options.min_control;
	}
}
void WorkingParam::ComputeRequiredBases(int NAA, int ss)
{
	// d: distance, fraction of errors;
	// e: number of errors;
	// g: length of the maximum gap;
	// m: word length;
	// n: sequence length;
	// alignment length = n - g + 1;
	// d = e / (n - g + 1);
	// e >= 1, so that, g <= n + 1 - 1/d
	// word count = (n - g - m + 1) - (e - 1)*m;
	//            = (n - g - m + 1) - (d*(n - g + 1) - 1)*m
	//            = (n - g + 1) - d*m*(n - g + 1)
	//            = (n - g + 1)*(1 - d*m)
	// minimum word count is reached when g == n + 1 - 1/d
	// so, minimum word count = 1/d - m.
	// if g == band_width: word count = (n - band + 1)*(1 - d*m);
	if (options.useDistance) {
		// int band = options.band_width + 1;
		int invd = int(1.0 / (options.distance_thd + 1E-9));
		// int k = len_eff < invd ? len_eff : invd;
		int ks = len_eff - ss + 1;
		int kn = len_eff - NAA + 1;
		int ks2 = invd - ss;
		int kn2 = invd - NAA;
		// int ks3 = int((len_eff - band + 1.0)*(1.0 - options.distance_thd * ss));
		// int kn3 = int((len_eff - band + 1.0)*(1.0 - options.distance_thd * NAA));
		//if( ks3 > ks2 ) ks2 = ks3;
		//if( kn3 > kn2 ) kn2 = kn3;
		required_aa1 = required_aas = (ks2 < ks ? ks2 : ks);
		required_aan = kn2 < kn ? kn2 : kn;
		if (required_aa1 <= 0) required_aa1 = required_aas = 1;
		if (required_aan <= 0) required_aan = 1;
		//required_aa1 = required_aas = required_aan = 0;
		return;
	}
	// (N-K)-K*(1-C)*N = C*K*N-(K-1)*N-K = (C*K-K+1)*N-K
	required_aa1 = (len_eff - ss) - int(ss * ceil((1.0 - aa1_cutoff) * len_eff));
	if (required_aa1 < 0) required_aa1 = 0;
	required_aas = required_aa1;
	required_aan = (len_eff - NAA) - int(NAA * ceil((1.0 - aa1_cutoff) * len_eff));
	//printf( "%i %i\n", required_aa1, required_aan );
	if (required_aan < 0) required_aan = 0;

	int aa1_old = int(aa1_cutoff * (double)len_eff) - ss + 1;
	int aas_old = int(aas_cutoff * (double)len_eff);
	int aan_old = int(aan_cutoff * (double)len_eff);

	double thd = option.cluster_thd;
	//double rest = (len_eff - ss) / double(len_eff * ss);
	double rest = (len_eff - NAA) / double(len_eff * NAA);
	double thd0 = 1.0 - rest;
	double fnew = 0;
	double fold = 1;
	if (thd > thd0) {
		fnew = (thd - thd0) / rest;
		fold = 1.0 - fnew;
	}
	//printf( "%g %g %g\n", thd, thd0, fnew );

	required_aa1 = (int)(fnew * required_aa1 + fold * aa1_old);
	required_aas = (int)(fnew * required_aas + fold * aas_old);
	required_aan = (int)(fnew * required_aan + fold * aan_old);
}


// when alignment coverage such as -aL is specified
// if a existing rep is too long, it won't be qulified 




size_t MemoryLimit( size_t mem_need )
{
	size_t mem_limit = (options.max_memory - mem_need) / sizeof(IndexCount);

	//printf( "Table limit with the given memory limit:\n" );
	if( options.max_memory == 0 ){
		mem_limit = options.max_entries;
		if( mem_limit > MAX_TABLE_SIZE ) mem_limit = MAX_TABLE_SIZE;
	}
	//printf( "Max number of representatives: %zu\n", mem_limit );
	//printf( "Max number of word counting entries: %zu\n\n", mem_limit );
	return mem_limit;
}

int calc_ann_list(int len, char *seqi, int NAA, int& aan_no, std::vector<int> & aan_list, std::vector<INTs> & aan_list_no, bool est) 
{
	// check_aan_list 
	const auto aan_no = len - NAA + 1;
	for (j=0; j<aan_no; j++) {
		aan_list[j] = 0;
		for (k=0, k1=NAA-1; k<NAA; k++, k1--) aan_list[j] += seqi[j+k] * NAAN_array[k1];
	}
	if( est ){
		// for the short word containing 'N', mask it to '-1'
		for (j=0; j<len; j++){
			if ( seqi[j] >= 4 ) {                      // here N is 4
				i0 = (j-NAA+1 > 0) ? j-NAA+1 : 0;
				i1 = j < aan_no ? j : aan_no - 1;
				for (i=i0; i<=i1; i++) aan_list[i]=-1;
			}
		}
	}

	std::sort(aan_list);

	for(j=0; j<aan_no; j++)
		aan_list_no[j]=1;

	for(j=aan_no-1; j; j--) {
		if (aan_list[j] == aan_list[j-1]) {
			aan_list_no[j-1] += aan_list_no[j];
			aan_list_no[j]=0;
		}
	}
	return OK_FUNC;
} // END calc_ann_list

//quick_sort_idx calling (a, idx, 0, no-1)
//sort a with another array idx
//so that idx rearranged
int quick_sort_idx (int *a, int *idx, int lo0, int hi0 ) {
  int lo = lo0;
  int hi = hi0;
  int mid;
  int tmp;

  if ( hi0 > lo0) {
    mid = a[ ( lo0 + hi0 ) / 2 ];

    while( lo <= hi ) {
      while( ( lo < hi0 ) && ( a[lo] < mid ) ) lo++;
      while( ( hi > lo0 ) && ( a[hi] > mid ) ) hi--;
      if( lo <= hi ) {
        tmp=a[lo];   a[lo]=a[hi];     a[hi]=tmp;
        tmp=idx[lo]; idx[lo]=idx[hi]; idx[hi]=tmp;
        lo++; hi--;
      }
    } // while

    if( lo0 < hi ) quick_sort_idx(a, idx, lo0, hi );
    if( lo < hi0 ) quick_sort_idx(a, idx, lo, hi0 );
  } // if ( hi0 > lo0)
  return 0;
} // quick_sort_idx


//decreasing can not use reverse of quick_sort_idx due to tie
//quick_sort_idxr calling (a, idx, 0, no-1)
//sort a with another array idx
//so that idx rearranged
int quick_sort_idxr (int *a, int *idx, int lo0, int hi0 ) {
  int lo = lo0;
  int hi = hi0;
  int mid;
  int tmp;

  if ( hi0 > lo0) {
    mid = a[ ( lo0 + hi0 ) / 2 ];

    while( lo <= hi ) {
      while( ( lo < hi0 ) && ( a[lo] > mid ) ) lo++;
      while( ( hi > lo0 ) && ( a[hi] < mid ) ) hi--;
      if( lo <= hi ) {
        tmp=a[lo];   a[lo]=a[hi];     a[hi]=tmp;
        tmp=idx[lo]; idx[lo]=idx[hi]; idx[hi]=tmp;
        lo++; hi--;
      }
    } // while

    if( lo0 < hi ) quick_sort_idxr(a, idx, lo0, hi );
    if( lo < hi0 ) quick_sort_idxr(a, idx, lo, hi0 );
  } // if ( hi0 > lo0)
  return 0;
} // quick_sort_idxr

/////////////////////////// END ALL ////////////////////////

int naa_stat_start_percent = 40;
int naa_stat[5][61][4] = {

	// cover 0.99
	{
		// N=5   N=4   N=3   N=2
		{  0,    0,    0,    7,  },  // 40%
		{  0,    0,    0,    8,  },  // 41%
		{  0,    0,    0,    9,  },  // 42%
		{  0,    0,    0,    9,  },  // 43%
		{  0,    0,    1,   10,  },  // 44%
		{  0,    0,    1,   11,  },  // 45%
		{  0,    0,    1,   12,  },  // 46%
		{  0,    0,    2,   13,  },  // 47%
		{  0,    0,    2,   14,  },  // 48%
		{  0,    0,    4,   16,  },  // 49%
		{  0,    0,    4,   16,  },  // 50%
		{  0,    0,    5,   17,  },  // 51%
		{  0,    0,    5,   18,  },  // 52%
		{  0,    0,    7,   20,  },  // 53%
		{  0,    1,    7,   21,  },  // 54%
		{  0,    1,    7,   21,  },  // 55%
		{  0,    2,    8,   23,  },  // 56%
		{  0,    2,    8,   25,  },  // 57%
		{  0,    2,   10,   25,  },  // 58%
		{  0,    3,   10,   26,  },  // 59%
		{  0,    4,   13,   28,  },  // 60%
		{  0,    5,   13,   30,  },  // 61%
		{  0,    5,   14,   30,  },  // 62%
		{  1,    6,   15,   33,  },  // 63%
		{  2,    7,   17,   34,  },  // 64%
		{  2,    7,   17,   35,  },  // 65%
		{  2,    9,   20,   37,  },  // 66%
		{  4,   10,   20,   37,  },  // 67%
		{  4,   11,   22,   40,  },  // 68%
		{  5,   12,   24,   41,  },  // 69%
		{  5,   12,   25,   42,  },  // 70%
		{  6,   16,   27,   43,  },  // 71%
		{  8,   16,   27,   45,  },  // 72%
		{  9,   17,   29,   47,  },  // 73%
		{ 10,   18,   31,   47,  },  // 74%
		{ 10,   20,   32,   50,  },  // 75%
		{ 12,   20,   32,   51,  },  // 76%
		{ 14,   22,   36,   54,  },  // 77%
		{ 15,   24,   37,   55,  },  // 78%
		{ 17,   26,   41,   58,  },  // 79%
		{ 18,   29,   41,   59,  },  // 80%
		{ 20,   30,   45,   60,  },  // 81%
		{ 24,   35,   48,   62,  },  // 82%
		{ 26,   36,   48,   64,  },  // 83%
		{ 27,   38,   51,   65,  },  // 84%
		{ 31,   43,   54,   68,  },  // 85%
		{ 35,   43,   55,   70,  },  // 86%
		{ 36,   48,   60,   71,  },  // 87%
		{ 36,   50,   61,   73,  },  // 88%
		{ 40,   50,   61,   75,  },  // 89%
		{ 45,   54,   65,   75,  },  // 90%
		{ 52,   60,   70,   79,  },  // 91%
		{ 53,   62,   71,   81,  },  // 92%
		{ 57,   66,   75,   84,  },  // 93%
		{ 57,   66,   76,   85,  },  // 94%
		{ 64,   71,   78,   85,  },  // 95%
		{ 70,   75,   82,   89,  },  // 96%
		{ 77,   81,   86,   92,  },  // 97%
		{ 82,   86,   90,   94,  },  // 98%
		{ 83,   87,   91,   95,  },  // 99%
		{ 91,   93,   95,   97,  },  // 100%
	},
	// cover 0.95
	{
		// N=5   N=4   N=3   N=2
		{  0,    0,    1,    9,  },  // 40%
		{  0,    0,    2,   10,  },  // 41%
		{  0,    0,    2,   11,  },  // 42%
		{  0,    0,    3,   12,  },  // 43%
		{  0,    0,    3,   12,  },  // 44%
		{  0,    0,    4,   14,  },  // 45%
		{  0,    0,    4,   14,  },  // 46%
		{  0,    1,    5,   16,  },  // 47%
		{  0,    1,    6,   17,  },  // 48%
		{  0,    2,    7,   19,  },  // 49%
		{  0,    2,    8,   19,  },  // 50%
		{  0,    2,    8,   20,  },  // 51%
		{  0,    2,    9,   21,  },  // 52%
		{  0,    4,   10,   23,  },  // 53%
		{  1,    4,   11,   24,  },  // 54%
		{  1,    4,   11,   24,  },  // 55%
		{  1,    5,   13,   26,  },  // 56%
		{  2,    5,   13,   27,  },  // 57%
		{  2,    6,   15,   29,  },  // 58%
		{  2,    7,   15,   30,  },  // 59%
		{  3,    8,   16,   31,  },  // 60%
		{  4,    8,   18,   32,  },  // 61%
		{  4,    9,   18,   33,  },  // 62%
		{  5,   11,   20,   36,  },  // 63%
		{  6,   12,   22,   37,  },  // 64%
		{  6,   12,   22,   38,  },  // 65%
		{  8,   14,   24,   40,  },  // 66%
		{  8,   15,   25,   41,  },  // 67%
		{ 10,   16,   27,   42,  },  // 68%
		{ 10,   18,   28,   45,  },  // 69%
		{ 11,   18,   29,   45,  },  // 70%
		{ 14,   21,   31,   47,  },  // 71%
		{ 14,   22,   32,   48,  },  // 72%
		{ 14,   22,   33,   50,  },  // 73%
		{ 17,   24,   36,   52,  },  // 74%
		{ 17,   25,   36,   52,  },  // 75%
		{ 18,   27,   39,   54,  },  // 76%
		{ 20,   29,   41,   56,  },  // 77%
		{ 21,   31,   42,   58,  },  // 78%
		{ 21,   31,   46,   60,  },  // 79%
		{ 27,   35,   46,   60,  },  // 80%
		{ 28,   37,   50,   63,  },  // 81%
		{ 31,   38,   50,   64,  },  // 82%
		{ 34,   43,   53,   66,  },  // 83%
		{ 36,   45,   54,   67,  },  // 84%
		{ 41,   50,   60,   70,  },  // 85%
		{ 43,   51,   60,   71,  },  // 86%
		{ 45,   54,   63,   74,  },  // 87%
		{ 48,   55,   64,   75,  },  // 88%
		{ 54,   60,   68,   78,  },  // 89%
		{ 55,   62,   71,   80,  },  // 90%
		{ 56,   63,   71,   80,  },  // 91%
		{ 64,   70,   76,   84,  },  // 92%
		{ 69,   74,   80,   86,  },  // 93%
		{ 73,   78,   83,   88,  },  // 94%
		{ 74,   78,   84,   89,  },  // 95%
		{ 80,   84,   87,   91,  },  // 96%
		{ 83,   86,   90,   93,  },  // 97%
		{ 86,   89,   92,   95,  },  // 98%
		{ 91,   93,   95,   97,  },  // 99%
		{ 92,   93,   95,   97,  },  // 100%
	},
	// cover 0.9
	{
		// N=5   N=4   N=3   N=2
		{  0,    0,    2,   11,  },  // 40%
		{  0,    0,    3,   12,  },  // 41%
		{  0,    0,    3,   12,  },  // 42%
		{  0,    1,    4,   13,  },  // 43%
		{  0,    1,    5,   14,  },  // 44%
		{  0,    1,    5,   15,  },  // 45%
		{  0,    1,    6,   16,  },  // 46%
		{  0,    2,    7,   18,  },  // 47%
		{  0,    2,    7,   18,  },  // 48%
		{  0,    3,    9,   20,  },  // 49%
		{  1,    4,    9,   20,  },  // 50%
		{  1,    4,   10,   21,  },  // 51%
		{  1,    4,   11,   23,  },  // 52%
		{  2,    5,   12,   24,  },  // 53%
		{  2,    5,   12,   25,  },  // 54%
		{  2,    6,   13,   26,  },  // 55%
		{  3,    7,   14,   28,  },  // 56%
		{  3,    7,   15,   28,  },  // 57%
		{  4,    8,   16,   30,  },  // 58%
		{  5,    9,   17,   31,  },  // 59%
		{  5,   10,   18,   32,  },  // 60%
		{  6,   11,   20,   35,  },  // 61%
		{  6,   11,   20,   35,  },  // 62%
		{  7,   13,   22,   38,  },  // 63%
		{  8,   14,   23,   39,  },  // 64%
		{  8,   15,   24,   39,  },  // 65%
		{ 10,   16,   26,   42,  },  // 66%
		{ 10,   17,   27,   42,  },  // 67%
		{ 12,   19,   29,   44,  },  // 68%
		{ 13,   20,   30,   46,  },  // 69%
		{ 13,   21,   31,   47,  },  // 70%
		{ 16,   23,   33,   48,  },  // 71%
		{ 18,   25,   34,   50,  },  // 72%
		{ 18,   26,   36,   51,  },  // 73%
		{ 19,   28,   38,   53,  },  // 74%
		{ 20,   29,   38,   53,  },  // 75%
		{ 23,   30,   41,   56,  },  // 76%
		{ 24,   33,   43,   57,  },  // 77%
		{ 26,   34,   45,   59,  },  // 78%
		{ 28,   37,   48,   61,  },  // 79%
		{ 30,   37,   48,   62,  },  // 80%
		{ 33,   42,   52,   64,  },  // 81%
		{ 35,   43,   53,   65,  },  // 82%
		{ 38,   47,   56,   68,  },  // 83%
		{ 40,   47,   56,   68,  },  // 84%
		{ 44,   53,   61,   71,  },  // 85%
		{ 45,   53,   62,   73,  },  // 86%
		{ 50,   58,   66,   75,  },  // 87%
		{ 51,   58,   66,   76,  },  // 88%
		{ 57,   63,   71,   79,  },  // 89%
		{ 60,   66,   72,   81,  },  // 90%
		{ 62,   68,   75,   83,  },  // 91%
		{ 70,   74,   80,   85,  },  // 92%
		{ 74,   78,   82,   88,  },  // 93%
		{ 85,   87,   90,   92,  },  // 94%
		{ 86,   88,   90,   92,  },  // 95%
		{ 87,   89,   91,   93,  },  // 96%
		{ 87,   89,   92,   94,  },  // 97%
		{ 89,   91,   93,   96,  },  // 98%
		{ 93,   94,   96,   97,  },  // 99%
		{ 94,   95,   97,   98,  },  // 100%
	},
	// cover 0.8
	{
		// N=5   N=4   N=3   N=2
		{  0,    1,    4,   13,  },  // 40%
		{  0,    1,    5,   13,  },  // 41%
		{  0,    1,    5,   14,  },  // 42%
		{  0,    2,    6,   15,  },  // 43%
		{  0,    2,    6,   16,  },  // 44%
		{  0,    2,    7,   17,  },  // 45%
		{  1,    3,    8,   18,  },  // 46%
		{  1,    4,    9,   20,  },  // 47%
		{  1,    4,    9,   20,  },  // 48%
		{  2,    5,   11,   22,  },  // 49%
		{  2,    5,   11,   22,  },  // 50%
		{  2,    6,   12,   24,  },  // 51%
		{  3,    6,   13,   25,  },  // 52%
		{  3,    7,   14,   26,  },  // 53%
		{  4,    8,   14,   27,  },  // 54%
		{  4,    8,   15,   28,  },  // 55%
		{  5,    9,   17,   30,  },  // 56%
		{  5,    9,   17,   30,  },  // 57%
		{  6,   11,   19,   32,  },  // 58%
		{  7,   12,   20,   34,  },  // 59%
		{  8,   12,   20,   34,  },  // 60%
		{  9,   14,   22,   37,  },  // 61%
		{  9,   14,   23,   37,  },  // 62%
		{ 10,   16,   25,   39,  },  // 63%
		{ 11,   17,   26,   41,  },  // 64%
		{ 12,   18,   27,   41,  },  // 65%
		{ 13,   20,   28,   43,  },  // 66%
		{ 14,   21,   30,   45,  },  // 67%
		{ 15,   22,   31,   46,  },  // 68%
		{ 17,   24,   33,   48,  },  // 69%
		{ 17,   24,   34,   48,  },  // 70%
		{ 19,   26,   36,   50,  },  // 71%
		{ 20,   27,   37,   51,  },  // 72%
		{ 21,   29,   39,   53,  },  // 73%
		{ 23,   31,   41,   55,  },  // 74%
		{ 23,   31,   41,   55,  },  // 75%
		{ 26,   34,   44,   58,  },  // 76%
		{ 28,   36,   46,   59,  },  // 77%
		{ 29,   37,   47,   60,  },  // 78%
		{ 34,   41,   50,   62,  },  // 79%
		{ 34,   42,   51,   63,  },  // 80%
		{ 38,   45,   55,   66,  },  // 81%
		{ 39,   46,   55,   67,  },  // 82%
		{ 44,   51,   60,   70,  },  // 83%
		{ 44,   51,   60,   70,  },  // 84%
		{ 49,   56,   64,   73,  },  // 85%
		{ 50,   57,   64,   74,  },  // 86%
		{ 57,   63,   69,   77,  },  // 87%
		{ 58,   64,   70,   78,  },  // 88%
		{ 68,   71,   76,   82,  },  // 89%
		{ 68,   72,   77,   83,  },  // 90%
		{ 75,   79,   81,   85,  },  // 91%
		{ 86,   87,   89,   90,  },  // 92%
		{ 88,   89,   90,   92,  },  // 93%
		{ 90,   91,   92,   93,  },  // 94%
		{ 91,   92,   93,   94,  },  // 95%
		{ 92,   94,   94,   95,  },  // 96%
		{ 93,   94,   95,   96,  },  // 97%
		{ 94,   95,   95,   96,  },  // 98%
		{ 94,   95,   96,   98,  },  // 99%
		{ 95,   96,   97,   98,  },  // 100%
	},
	// cover 0.6
	{
		// N=5   N=4   N=3   N=2
		{  1,    2,    6,   15,  },  // 40%
		{  1,    3,    7,   16,  },  // 41%
		{  1,    3,    8,   17,  },  // 42%
		{  2,    4,    9,   18,  },  // 43%
		{  2,    4,    9,   19,  },  // 44%
		{  2,    5,   10,   20,  },  // 45%
		{  3,    5,   10,   21,  },  // 46%
		{  3,    6,   12,   22,  },  // 47%
		{  3,    6,   12,   23,  },  // 48%
		{  4,    8,   14,   25,  },  // 49%
		{  4,    8,   14,   25,  },  // 50%
		{  5,    8,   15,   26,  },  // 51%
		{  5,    9,   16,   27,  },  // 52%
		{  6,   10,   17,   29,  },  // 53%
		{  6,   11,   18,   30,  },  // 54%
		{  7,   11,   18,   31,  },  // 55%
		{  8,   12,   20,   32,  },  // 56%
		{  8,   13,   20,   33,  },  // 57%
		{ 10,   14,   22,   35,  },  // 58%
		{ 10,   15,   23,   37,  },  // 59%
		{ 11,   16,   24,   37,  },  // 60%
		{ 12,   18,   26,   39,  },  // 61%
		{ 13,   18,   26,   40,  },  // 62%
		{ 14,   20,   28,   42,  },  // 63%
		{ 16,   22,   30,   43,  },  // 64%
		{ 16,   22,   31,   44,  },  // 65%
		{ 17,   23,   32,   45,  },  // 66%
		{ 18,   25,   33,   47,  },  // 67%
		{ 19,   26,   35,   48,  },  // 68%
		{ 21,   27,   36,   50,  },  // 69%
		{ 22,   29,   37,   51,  },  // 70%
		{ 24,   30,   39,   52,  },  // 71%
		{ 25,   32,   41,   53,  },  // 72%
		{ 26,   33,   42,   55,  },  // 73%
		{ 29,   35,   44,   57,  },  // 74%
		{ 29,   36,   45,   57,  },  // 75%
		{ 32,   39,   48,   60,  },  // 76%
		{ 34,   41,   50,   61,  },  // 77%
		{ 36,   43,   51,   62,  },  // 78%
		{ 40,   46,   54,   65,  },  // 79%
		{ 40,   46,   54,   65,  },  // 80%
		{ 46,   52,   59,   68,  },  // 81%
		{ 46,   52,   60,   69,  },  // 82%
		{ 53,   59,   65,   73,  },  // 83%
		{ 54,   60,   66,   73,  },  // 84%
		{ 63,   67,   73,   78,  },  // 85%
		{ 68,   71,   75,   79,  },  // 86%
		{ 78,   80,   82,   85,  },  // 87%
		{ 79,   81,   83,   85,  },  // 88%
		{ 83,   85,   86,   87,  },  // 89%
		{ 85,   86,   87,   89,  },  // 90%
		{ 86,   88,   89,   90,  },  // 91%
		{ 88,   89,   90,   91,  },  // 92%
		{ 90,   90,   91,   92,  },  // 93%
		{ 91,   92,   92,   93,  },  // 94%
		{ 92,   93,   94,   94,  },  // 95%
		{ 94,   94,   95,   95,  },  // 96%
		{ 95,   95,   96,   96,  },  // 97%
		{ 95,   96,   97,   97,  },  // 98%
		{ 96,   96,   97,   98,  },  // 99%
		{ 97,   98,   98,   99,  },  // 100%
	},
};

