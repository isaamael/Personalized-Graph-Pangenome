#include <Rcpp.h>
#include <algorithm>
#include <cmath>
using namespace Rcpp;
// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::plugins(openmp)]]
// Lossless run-boundary representation; prefix integrates every evaluated window.
inline double gd(int q,int r,const IntegerVector &p,const IntegerVector &ii,const NumericVector &xx,const NumericVector &pre){
 int a=p[q],ae=p[q+1],b=p[r],be=p[r+1],sa=0,sb=0,pos=0,n=pre.size()-1;double d=0;
 while(a<ae||b<be){int nx=std::min(a<ae?ii[a]:n,b<be?ii[b]:n);d+=(pre[nx]-pre[pos])*std::abs(sa-sb);pos=nx;if(a<ae&&ii[a]==nx){sa+=int(xx[a]);a++;}if(b<be&&ii[b]==nx){sb+=int(xx[b]);b++;}}
 return d+(pre[n]-pre[pos])*std::abs(sa-sb);
}
// [[Rcpp::export]]
NumericVector pair_distances(S4 X,IntegerVector a,IntegerVector b,NumericVector prefix,int threads=4){
 IntegerVector p=X.slot("p"),ii=X.slot("i");NumericVector xx=X.slot("x"),out(a.size());
 #pragma omp parallel for num_threads(threads) schedule(static)
 for(int k=0;k<a.size();k++)out[k]=gd(a[k]-1,b[k]-1,p,ii,xx,prefix);return out;
}
// Pairwise SGD minimizes unweighted raw metric stress; sampled edges approximate
// the full all-observation objective. No route, generation, or GEBV enters fitting.
// [[Rcpp::export]]
NumericMatrix stress_sgd(NumericMatrix initial,IntegerVector a,IntegerVector b,NumericVector d,int epochs=100,double first=.5,double last=.0001){
 NumericMatrix z=clone(initial);int m=a.size();
 for(int e=0;e<epochs;e++){double rate=first*std::pow(last/first,double(e)/std::max(1,epochs-1));int start=int((static_cast<long long>(e)*104729)%m);
  for(int t=0;t<m;t++){int k=(t+start)%m,i=a[k]-1,j=b[k]-1;if(i==j)continue;double dx=z(i,0)-z(j,0),dy=z(i,1)-z(j,1),r=std::sqrt(dx*dx+dy*dy);if(r<1e-12)continue;double f=.5*rate*(r-d[k])/r;z(i,0)-=f*dx;z(i,1)-=f*dy;z(j,0)+=f*dx;z(j,1)+=f*dy;}
  if(e%20==0)Rcpp::checkUserInterrupt();
 }return z;
}
// Each query independently minimizes squared distance error to ALL fixed F2s.
// Three deterministic starts; gradient majorization with line search.
// [[Rcpp::export]]
List project_fixed(S4 X,IntegerVector queries,IntegerVector anchors,NumericMatrix az,NumericMatrix initial,NumericVector prefix,int threads=4){
 IntegerVector p=X.slot("p"),ii=X.slot("i");NumericVector xx=X.slot("x");int nq=queries.size(),na=anchors.size();NumericMatrix z(nq,2);NumericVector err(nq);
 #pragma omp parallel for num_threads(threads) schedule(dynamic,8)
 for(int q=0;q<nq;q++){
  std::vector<double> ds(na);double mean=0;for(int j=0;j<na;j++){ds[j]=gd(queries[q]-1,anchors[j]-1,p,ii,xx,prefix);mean+=ds[j]/na;}
  auto loss=[&](double x,double y){double v=0;for(int j=0;j<na;j++){double dx=x-az(j,0),dy=y-az(j,1),r=std::sqrt(dx*dx+dy*dy)-ds[j];v+=r*r;}return v/na;};
  double best=1e100,bx=0,by=0;
  for(int s=0;s<3;s++){
   double x=s==0?initial(q,0):((s==1?.7:-.7)*mean),y=s==0?initial(q,1):.4*mean,prev=loss(x,y);
   for(int it=0;it<200;it++){double gx=0,gy=0;for(int j=0;j<na;j++){double dx=x-az(j,0),dy=y-az(j,1),r=std::max(1e-12,std::sqrt(dx*dx+dy*dy)),f=(r-ds[j])/r;gx+=f*dx/na;gy+=f*dy/na;}
    double step=1,cur=prev,nx=x,ny=y;for(int ls=0;ls<15;ls++){nx=x-step*gx;ny=y-step*gy;cur=loss(nx,ny);if(cur<=prev)break;step*=.5;}
    double change=prev-cur;x=nx;y=ny;prev=cur;if(change>=0&&change<1e-10)break;
   }
   if(prev<best){best=prev;bx=x;by=y;}
  }z(q,0)=bx;z(q,1)=by;err[q]=best;
 }return List::create(_["xy"]=z,_["mean_squared_distance_error"]=err);
}
