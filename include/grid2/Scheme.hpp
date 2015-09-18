#ifndef SCHEME
#define SCHEME

class Var; 

template <typename T> class Scheme {
public:
  vector<int_8> ind; 
  vector<T> val; 
  T c; 
  int_2 type;  // 0: cell center, 1: vertex, 2: face center; 

  Scheme(): type(0), c(0.0) {}
  Scheme(int_2 t): type(t), c(0.0) {}

  void resize(int n) {
    ind.resize(n); val.resize(n); 
  }

  int size() { 
    return ind.size(); 
  }

  bool isDuplicate(int_8 i) {
    return (std::find(ind.begin(), ind.end(), i) != ind.end()); 
  }

  // int_8 varSize() {
  //   if (type == 0) return listCell.size(); 
  //   else if (type == 1) return listVertex.size(); 
  //   else if (type == 2) return listFace.size(); 
  //   else {
  //     cout << "Scheme:varSize(): type is not understood! "<< type << endl; 
  //     return 0; 
  //   }
  // }    

  void push_pair(int_8 i, T v) {
    auto it = std::find(ind.begin(), ind.end(), i); 
    if (it == ind.end()) {
      ind.emplace_back(i); 
      val.emplace_back(v); 
    } else {
      auto i = it - ind.begin(); 
      val[i] += v; 
    }
  }
  void push_constant(T d) {
    c += d; 
  }

  T eval(shared_ptr<Var > phi) {
    return eval(phi->data); 
  }

  T eval(VecX<double> phi) {
    // if (phi.size() != varSize()) {
    //   cout << "Scheme:eval(): Size and location pointer do not match"<<endl; 
    //   exit(1); 
    // }
    T s=c; 
    for (auto i = 0; i < size(); ++i) {
      s += phi[ind[i]]*val[i]; 
    }
    return s ; 
  }


  Scheme<T>& operator/=(T a) { 
    if (a == 0) {cout <<  "Scheme:operator/: div by zero!"<< endl; exit(1);}
    for (auto i =0; i< size(); ++i) val[i] /= a; 
    return *this;
  } 

  Scheme<T>& operator*=(T a) { 
    if (a == 0) {cout <<  "Scheme:operator*: mult by zero!"<< endl; exit(1);}
    for (auto i =0; i< size(); ++i) val[i] *= a; 
    return *this;
  } 

  Scheme<T>& operator=(T const &a) { 
    ind.assign(a.ind.begin(), a.ind.end()); 
    val.assign(a.val.begin(), a.val.end()); 
    c = a.c; 
    return *this;
  } 

  Scheme<T>& operator+=(Scheme<T> const &a) {
    for (auto i =0; i < a.ind.size(); ++i) {
      push_pair(a.ind[i], a.val[i]); 
    }
    push_constant(a.c); 
    return *this;
  } 



};


#endif
