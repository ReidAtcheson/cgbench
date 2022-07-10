#include <vector>
#include <tuple>
#include <cmath>
#include <iostream>

using size_t = std::size_t;
//Makes a CSR format sparse matrix out of a 3D laplacian
std::tuple<std::vector<size_t>,std::vector<size_t>,std::vector<double>> make_matrix(size_t mx,size_t my,size_t mz){
  std::vector<size_t> cids;
  std::vector<size_t> offs;
  std::vector<double> vals;
  auto id = [&](size_t ix,size_t iy,size_t iz){
    return iz + mz*(iy + my*ix);
  };
  size_t off=0;
  offs.push_back(off);

  for(size_t ix=0;ix<mx;ix++){
    for(size_t iy=0;iy<my;iy++){
      for(size_t iz=0;iz<mz;iz++){
        auto add_nonzero = [&](size_t ix,size_t iy,size_t iz,double val){
          vals.push_back(val);
          cids.push_back(id(ix,iy,iz));
          off+=1;
        };
        add_nonzero(ix,iy,iz,6.0);
        if(ix>0){
            add_nonzero(ix-1,iy,iz,-1.0);
        }
        if(ix<mx-1){
            add_nonzero(ix+1,iy,iz,-1.0);
        }
        if(iy>0){
            add_nonzero(ix,iy-1,iz,-1.0);
        }
        if(iy<my-1){
            add_nonzero(ix,iy+1,iz,-1.0);
        }
        if(iz>0){
            add_nonzero(ix,iy,iz-1,-1.0);
        }
        if(iz<mz-1){
            add_nonzero(ix,iy,iz+1,-1.0);
        }
        offs.push_back(off);
      }
    }
  }
  return std::make_tuple(cids,offs,vals);
}


void csr_eval(const std::vector<size_t>& cids,const std::vector<size_t>& offs,const std::vector<double>& vals,const std::vector<double>& x,std::vector<double>& y){
  for(size_t r=0;r<y.size();r++){
    size_t beg = offs[r];
    size_t end = offs[r+1];
    y[r]=0.0;
    for(size_t k=beg;k<end;k++){
      size_t ri=cids[k];
      double val=vals[k];
      y[r] += val * x[ri];
    }
  }
}

double dot(const std::vector<double>& x,const std::vector<double>& y){
  double out=0.0;
  for(size_t i=0;i<x.size();i++){
    out+=x[i]*y[i];
  }
  return out;
}

void cgsolve(const std::vector<size_t>& cids,
    const std::vector<size_t>& offs,
    const std::vector<double>& vals,
    const std::vector<double>& b,
    std::vector<double>& r,
    std::vector<double>& p,
    std::vector<double>& ap,
    std::vector<double>& x){

  double tol = 1e-8;
  //Set r = b - Ax
  csr_eval(cids,offs,vals,x,r);
  for(size_t i=0;i<r.size();i++){
    r[i] = b[i] - r[i];
  }
  //Set p = r
  for(size_t i=0;i<p.size();i++){
    p[i]=r[i];
  }
  double rho = dot(r,r);
  while(std::sqrt(rho) > tol){
    csr_eval(cids,offs,vals,p,ap);
    double rdotr = rho;
    double pdotap = dot(p,ap);
    double alpha = rdotr/pdotap;
    for(size_t i=0;i<p.size();i++){
      x[i]+=alpha*p[i];
    }
    for(size_t i=0;i<p.size();i++){
      r[i]-=alpha*ap[i];
    }
    double rpdotrp = dot(r,r);
    double beta = rpdotrp/rdotr;
    for(size_t i=0;i<p.size();i++){
      p[i]=r[i]+beta*p[i];
    }
    rho=dot(r,r);
    std::cout<<"res = "<<std::sqrt(rho)<<"\n";
  }

}



int main(int argc,char** argv){
  if(argc==4){
    size_t mx = std::stoi(argv[1]);
    size_t my = std::stoi(argv[2]);
    size_t mz = std::stoi(argv[3]);
    std::tuple<std::vector<size_t>,std::vector<size_t>,std::vector<double>> mat = make_matrix(mx,my,mz);
    std::vector<double> x(mx*my*mz,0.0);
    std::vector<double> b(mx*my*mz,1.0);
    std::vector<double> r(mx*my*mz,0.0);
    std::vector<double> p(mx*my*mz,0.0);
    std::vector<double> ap(mx*my*mz,0.0);
    cgsolve(std::get<0>(mat),std::get<1>(mat),std::get<2>(mat),b,r,p,ap,x);
  }
  else{
    std::cout<<"Usage: "<<argv[0]<<" <mx> <my> <mz>\n";
    return 1;
  }
}


