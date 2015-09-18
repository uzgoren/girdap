#ifndef MATX
#define MATX

//-------------------------------------------
// Matrix storage scheme
//    Defines: basic matrix operations (not all)
//             Vec*Matrix addition
//-------------------------------------------



#include <VecX>

template <typename T> class MatX {
public: 
  //----------------------
  // Properties
  //----------------------
  vector<VecX<T> > data;
  int_8 rows, cols;  
  typedef typename vector<VecX<T> >::iterator rowiter;
  typedef typename vector<VecX<T> >::const_iterator crowiter;  
  
  typedef typename vector<T>::iterator it_vec; 
  typedef typename vector<T>::const_iterator it_cvec;
  typedef typename map<int_8, T>::iterator it_map;   
  typedef typename map<int_8, T>::const_iterator it_cmap;   

  //----------------------
  // Constructors
  //----------------------
  MatX() {rows=0; cols=0;}; 
  MatX(int_8 const &a) {
    resize(a); 
  }

  //----------------------
  // Methods
  //----------------------
  void reserve(int_8 const &a) {data.reserve(a);}
  void resize(int_8 const &a) {
    rows=a; cols=a; 
    if (data.capacity() == 0) data.reserve(a); 
    data.resize(a);
    int_8 row = 0; 
    for (rowiter it=data.begin(); it!=data.end(); ++it, ++row) {
      it->rank = cols; it->push_pair({row}, {0}); 
    }
  }
  void empty() {data.empty();}
  //Diagonal matrix assignment; b is the offset; 
  void insertDiag(VecX<T> const &a, int const &b=0) {
    if (a.rank != rows) {
      cout<< "Matx assignment by VecX : not same size! "<< a.rank << " " << rows <<endl; 
      exit(1); 
    }
    if (!a.px.empty()) {
      for (it_cmap it=a.px.begin(); it!=a.px.end(); ++it)
	if (b >= -it->first &&  b < a.rank-it->first)
	  data[it->first][it->first+b] = it->second;
    } else if (!a.data.empty()) {
      int_8 r = 0; 
      for (it_cvec it=a.data.begin(); it!=a.data.end(); ++it, ++r) {
	if (b >= (-r) && b < (a.rank - r)) 
	  data[r][r+b] = *it;
	}
    }
  }
  // Get diagonal elements; 
  VecX<T> getDiag() {
    VecX<T> a(rows); a.uncompress(); 
    for (int_8 r=0; r<rows; ++r) 
      a[r] = data[r][r];
    return a; 
  }

// Get upper diagonal elements; 
  VecX<T> getUpperDiag() {
    VecX<T> a(rows); a.uncompress(); 
    for (int_8 r=1; r<rows-1; ++r) 
      a[r] = data[r][r+1];
      a[0] =data[0][1];
     // a[rows] =data[rows-1][rows];
    return a; 
  }

// Get Lower diagonal elements; 
  VecX<T> getLowerDiag() {
    VecX<T> a(rows); a.uncompress(); 
    for (int_8 r=1; r<rows-1; ++r) 
      a[r] = data[r+1][r];
      a[0] =data[1][0];
      //a[rows] =data[rows][rows-1];
    return a; 
  }

  void info() {
    int_8 comp=0;  
    int_8 icnt=0; 
    int_8 bytes=0; 
    for (auto i = data.begin(); i !=data.end(); ++i) {
      icnt = icnt + (!(i->px.empty()) ? i->px.size() : i->data.size()); 
      bytes = bytes + (i->data.size())*sizeof(typename vector<double>::value_type) 
	+ (i->px.size())*sizeof(typename map<int_8,double>::value_type); 
      comp = comp + (i->px.empty() ? 0 : 1); 
    }
    cout<<endl<<"--------------------------------"<<endl; 
    cout<<"      Rows:\t"<<data.size()<<endl; 
    cout<<"    Length:\t"<<icnt <<endl; 
    cout<<"Total size:\t"<<bytes<<endl; 
    cout<<"Compressed:\t"<<comp<<endl; 
    cout<<"--------------------------------"<<endl; 
  }

  void compress() {
    for (auto i = data.begin(); i !=data.end(); ++i) i->compress();
  }

  void uncompress() {
    for (auto i = data.begin(); i !=data.end(); ++i) i->uncompress();
  }


  //----------------------
  // Operators
  //----------------------
  //01.Assignment! WITH MATRIX
  MatX<T>& operator=(MatX<T> const &a){
    this->data.assign(a.data.begin(),a.data.end());
    this->rows = a.rows; 
    this->cols = a.cols; 
    return *this;
  }

  //02. Assignment! with init_list
  MatX<T>& operator=(std::initializer_list<T> const &a) {
    int_8 row=0; 
    for (typename initializer_list<T>::const_iterator rit=a.begin(); rit!=a.end(); ++rit, ++row) {
      data[row][row] = *rit;
    }
    this->rows = row; this->cols = row; 
    return *this; 
  }

  //03. Assignment with double init_list: i.e. { {row1}, {row2}, ..., {rowN} }
  MatX<T>& operator=(std::initializer_list<initializer_list<T> > const &a){
    int_8 row=0; 
    for(typename initializer_list<initializer_list<T> >::const_iterator rit=a.begin(); rit!=a.end(); ++rit, ++row) {
      int_8 col=0; 
      for (typename initializer_list<T>::const_iterator cit=rit->begin(); cit!=rit->end(); ++cit, ++col) { 
	data[row][col] = *cit; 
      }
    }
    this->rows = row; this->cols = row; 
    return *this; 
  }

  VecX<T>& operator[] (int_8 const &index) {
    return data[index]; 
  };

  // O1: adds two matrix (expands rows if necessary - ISNEEDED?)
  MatX<T>& operator +=(MatX<T> const &a) {
    if (rows < a.rows) resize(a.rows);
    for (int_8 r = 0; r<rows; ++r) 
      if (r < a.data.size()) 
	data[r] = data[r] + a.data[r]; 
    return *this; 
  }

  MatX<T> operator +(const MatX<T>& a) {
    if (rows < a.rows) resize(a.rows);
    MatX<T> c(rows); 
    for (int_8 r = 0; r<rows; ++r) 
      //if (r < a.data.size()) 
  	c.data[r] = data[r] + a.data[r]; 
    return c; 
  }

  MatX<T> operator -(const MatX<T>& a) {
    if (rows < a.rows) resize(a.rows);
    MatX<T> c(rows); 
    for (int_8 r = 0; r<rows; ++r) 
      // if (r < a.data.size()) 
  	c.data[r] = data[r] - a.data[r]; 
    return c; 
  }

  MatX<T> operator -() {
    for (int_8 r=0; r < rows; ++r) 
      data[r] = -data[r]; 
    return (*this);
  }
  
  friend MatX<T> operator*( MatX<T> &a, T &b) {
    MatX<T> result(a.rows); result.compress(); 
    if (b == 0) return result; 
     for (int_8 r = 0; r < a.rows; ++r)
      result.data[r] = b*(a.data[r]); 
     return result;
  }

  friend MatX<T> operator*(T &b, MatX<T> &a) {return a*b;} 


  friend VecX<T> operator*(MatX<T>  &a, VecX<T>  &b) {
    VecX<T> result(a.rows); result.data.resize(a.rows); //uncompressed
    int_8 r=0; 
    //    cout << a << endl; 
    for (crowiter it=a.data.begin(); it!=a.data.end(); ++it, ++r) {
      //cout << "row (" << it->rank << "): " << *it << endl; 
      //for (r = 0; r < a.rows; ++r) {
      auto q = *it; //a[r]; //*(a.data.begin()+r); 
      result[r] = q*b;
    }
    return result; 
  }
  friend VecX<T> operator*(VecX<T> &b, MatX<T> const &a) {
    VecX<T> result(a.rows); 
    int_8 r=0; 
    for (crowiter it=a.data.begin(); it!=a.data.end(); ++it, ++r) {
      if (!it->px.empty()) {
	for (it_cmap cit=it->px.begin(); cit!=it->px.end(); ++cit) {
	  result[cit->first] += b[r]*cit->second; 
	}
      } else if(!it->data.empty()) {
	int_8 c=0; 
	for (it_cvec cit=it->data.begin(); cit!=it->data.end(); ++cit, ++c) {
	  result[c] += b[r]*(*cit); 
	}
      } 
    }
    return result; 
  }

  //----------------------
  // Output options
  //----------------------
  friend ostream &operator<<(ostream &out, const MatX<T> &a){
    int_8 row = 0; 
    for (crowiter it=a.data.begin(); it!=a.data.end(); ++it, ++row) 
      out << row << ": "<< *it << endl;
    return out;
  }

};



#endif
