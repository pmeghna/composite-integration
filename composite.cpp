/*
(c) M. Alamgir Hossain. Using this program I solved composite numerical integration. Which already done in our following paper using Mathematica code: 

Md. Shafiqul Islam and M. Alamgir Hossain (2009)- Numerical Integrations over an Arbitrary Quadrilateral Region, Applied Mathematics and Computations Volume 210, Issue 2, 15 April 2009, Pages 515-524.
*/

#include<iostream>
#include<stdio.h>
#include<cmath>
#include<vector>

double f(double, double);
double F(double,double);
using namespace std;
int main()
{
	vector<double>c(50),x(50),y(50),h(50),App(51),xi(7),w(7);
	vector<double>::iterator it;

	double sum;
	int i,j,n,k,l,m;

	cout<<"Enter the order of Gauss Legendre Quadrature rule:";
	cin>>n;
	if(n==2){
		double x2[2]={-0.577350269189626,0.577350269189626};
		double w2[2]={1.0,1.0};
		xi.assign (x2,x2+n);
		w.assign(w2,w2+n);
	}

	else if(n==3) {
		double x3[3]={-0.774596669241483,0.0,0.774596669241483};
		double w3[3]={5.0/9.0,8.0/9.0,5.0/9.0};
		xi.assign (x3,x3+n);
		w.assign(w3,w3+n);
	}

	else if(n==4){
		double x4[4]={-0.861136311594053,-0.339981043584856,0.339981043584856,0.861136311594053};
		double w4[4]={0.347854845137454,0.652145154862546,0.652145154862546,0.347854845137454};
		xi.assign (x4,x4+n);
		w.assign(w4,w4+n);
	}
	else if(n==5){
		double x5[5]={-0.906179845938664,-0.538469310105683,0.0,0.538469310105683,0.906179845938664};
		double w5[5]={0.2369268850561890,0.478628670499366,0.568888888888889,0.478628670499366,0.2369268850561890};
		xi.assign (x5,x5+n);
		w.assign(w5,w5+n);
	}
	else if(n==6){
		double x6[6]={-0.932469514203152,-0.661209386466265,-0.238619186083197,0.238619186083197,0.661209386466265,0.932469514203152};
		double w6[6]={0.171324492379170,0.360761573048139,0.467913934572691,0.467913934572691,0.360761573048139,0.171324492379170};
		xi.assign (x6,x6+n);
		w.assign(w6,w6+n);
	}
	else if(n==7){
		double x7[7]={-0.949107912342759,-0.741531185599394,-0.405845151377397,0.0,0.405845151377397,0.741531185599394,0.949107912342759};
		double w7[7]={0.129484966168870,0.279705391489277,0.381830050505119,0.417959183673469,0.381830050505119,0.279705391489277,0.129484966168870};
		xi.assign (x7,x7+n);
		w.assign(w7,w7+n);
	}
	else{
		cout<<"Order of Gauss Legendre Quadrature rule not support in this program. ";
		return 1;	
	}
	cout<<"\tCk"<<"\t\t\t"<<"Xk"<<"\t\t\t"<<"Yk"<<endl;

	for(i=0;i<n;i++){
		for(j=0;j<n;j++){
			k=i+(n-j-1)*n;
			c[k]=(2-xi[i]-xi[j])*w[i]*w[j]/4.0;
			x[k]=(-1+3*xi[i]-xi[j]-xi[i]*xi[j])/4.0;
			y[k]=(-1-xi[i]+3*xi[j]-xi[i]*xi[j])/4.0;
		}
	}
	for(k=0;k<n*n;k++){
		printf("%0.15lf\t%0.15lf\t%0.15lf\n",c[k],x[k],y[k]);
	}
	
	cout<<"Enter the value of n (between 1 and 50) then discretise T into 4 n^2 triangle:";
	cin>>m;
	for(l=1;l<=m;l++){
		App[l]=0.0;
		sum=0.0;
		for(k=0;k<n*n;k++){
			h[k]=0.0;
			for(i=0;i<=2*l-1;i++){
				for(j=0;j<=2*l-1-i;j++){
                        h[k]=h[k]+F((x[k]+2*(i-l)+1)/(2*l),(y[k]+2*(j-l)+1)/(2*l));
				}
			}
			for(i=0;i<=2*l-2;i++){
				for(j=0;j<=2*l-2-i;j++){
						h[k]=h[k]+F((-x[k]+2*(i-l)+1)/(2*l),(-y[k]+2*(j-l)+1)/(2*l));
				}
			}
			sum=sum+c[k]*h[k];
		}
		App[l]=sum/(4*l*l);
	}
	cout<<"4*n^2"<<"\t\t"<<"Result"<<endl;
	for(i=1;i<=m;i++){
		printf("4*%d\t%0.15lf\n",i*i,App[i]);
	}
	return 0;
}

/* Example 2 */
double f(double x, double y){
	double c1=65625.0/208.0,c2=328125.0/104.0,c3=239062.0/208.0,d1=1.0,d2=-125.0/4.0,d3=175.0/4.0;
	return (c1*pow(x,8)+c2*pow(y,9)+c3*pow(x,7)*pow(y,6))/(d1*pow(x,9)+d2*pow(y,7)+d3);
}

double F(double x,double y){
	return f(x,y)+f(-y,-x);

}
