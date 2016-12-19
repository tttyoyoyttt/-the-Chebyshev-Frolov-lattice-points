#include <stdio.h>
#include <math.h>


//This is the implementation of Algorithm 2 in the paper (arxiv:1612.05342)
//"Enumeration of the Chebyshev-Frolov lattice points in axis-parallel boxes",
//which counts the number of the Chebyshev-Frolov lattice points in [-1/2,1/2]^d for the scaling parameter N.
//This is made by Kosuke Suzuki and Takehito Yoshiki.
//Please read that article about the detailed information.

double Max(double a,double b){
	return a>b?a:b;
}
//The home made Maximum function: Max(x,y)


double Min(double a,double b){
	return a<b?a:b;
}
//The home made Minimum function: Min(x,y)


#define nmax 10
#define dmax 600
//Through this implementation,
//d denotes the dimension of the unit cube and n denotes log2( d ).
//We set nmax=10 such that n < nmax.
//We define dmax such that d=2^n< 2^(nmax-1)+1 <= dmax.

//We observed that, for a fixed nmax, the execution time varies depending on dmax.
//Thus we choose dmax such that the execution time is smaller heuristically.

//INPUT :
//arg[1]=n, arg[2]=log2( the scaling parameter N )
//OUTPUT:
//printf("dim=%d, log2(N)=%d, nodes=%lld\n")
//dim      = the dimension of the unit cube
//log2(N)  = log2 of the scaling parameter
//nodes    = the number of the Chebyshev-Frolov lattice points
//We count the cardinality of the nodes for the lattice points.
//We do not store the nodes.


double a[nmax][dmax];
double b[nmax][dmax];
double c[nmax][dmax];
//alpha[L][a] = (a[L][(a-1)*2^L+1],...,a[L][a*2^L]),
//beta [L][a] = (b[L][(a-1)*2^L+1],...,b[L][a*2^L]),
//gamma[L][a] = (c[L][(a-1)*2^L+1],...,c[L][a*2^L]).


int k[dmax];
//This array represents the integer point with d components (k[1],...,k[d]) for each point of the Chebyshev-Frolov lattice.


double	Cos[nmax][dmax];
double	Cosinv[nmax][dmax];
//D_n=[Cos[n][1],...,Cos[n][d]] and (D_n)^(-1)=[Cosinv[n][1],...,Cosinv[n][d]].
//We calculate (D_n)^(-1) in order to reduce the execution time.


int decomp[dmax][2];
//i is decomposed as i = 2^(decomp[i][0]) * decomp[i][1] where decomp[i][1] is the odd number.


long long int ptcount;
//ptcount = the number of the Chebyshev-Frolov lattice points.


int sigma(int n, int k){
	int d=1<<(n-1);
	if ((1<=k)&&(k<=d)) return sigma(n-1,k);
	else if ((d+1<=k)&&(k<=2*d)) return 2*d+1-sigma(n-1,k-d);
	else if ((n==0)&&(k==1)) return 1;
}
//This is the definition of sigma(n,k).


int Enum(int i,int n){
	//Enumerating the i-th coordinate k[i]
	for(k[i] = ceil(b[0][i]);k[i] <= floor(c[0][i]);k[i]++){
			if( i != (1<<n) ){
			int r=decomp[i][0],p=decomp[i][1]; //decompose i as i = 2^r*p
			a[0][i]=k[i];
			int j;
			for( j = 1; j<= r; j++){
				//Updating a[L][a]
				int h;
				for(h=1; h<=(1<<(j-1)); h++ ){
					int preind1=1<<(j-1);
					int preind2=(p*(1<<(r-j+1))-2)*preind1+h;
					int preind3=preind1+preind2;
					double a1=a[j-1][preind2];
					double a2=a[j-1][preind3]*Cos[j-1][h];
					a[j][preind2]=a1 + a2;
					a[j][preind3]=a1 - a2;
				}				
			}
			int preind1=(1<<r);
			int h;
			for(h=1; h<=preind1; h++ ){
				int preind2=p*preind1+h;
				int preind3=preind2-preind1;
				b[r][preind2]=Max(b[r+1][preind3]-a[r][preind3],-c[r+1][preind2]+a[r][preind3])*Cosinv[r][h];
				c[r][preind2]=Min(c[r+1][preind3]-a[r][preind3],-b[r+1][preind2]+a[r][preind3])*Cosinv[r][h];
			}
			for(j = r-1; j>=0; j--){
				//Updating b[L][a] and c[L][a]
				int preind1=(1<<j);
				int h;
				for (h=1; h<=(1<<j); h++ ){
					int preind2=p*(1<<r)+h;
					int preind3=preind1+preind2;
					b[j][preind2]=(b[j+1][preind2]+b[j+1][preind3])*0.5;
					c[j][preind2]=(c[j+1][preind2]+c[j+1][preind3])*0.5;
				}
			}
			Enum(i + 1,n);
		}else{
			//That is, if i = 2^n
			a[0][i]=k[i];
			int preind1=(1<<n);
			int j;
			for( j = 1; j <= n; j++ ){
				//Updating a[n][a]
				int preind2=(1<<(j-1));
				int h;
				for (h=1; h<=(1<<(j-1)); h++) {
					int preind3=preind1-preind2+h;
					int preind4=preind3-preind2;
					double a1=a[j-1][preind4];
					double a2=a[j-1][preind3]*Cos[j-1][h];
					a[j][preind4]=a1 + a2;
					a[j][preind3]=a1 - a2;
				}				
			}
		ptcount++; //Counting the number of the Chebyshev-Frolov lattice points
		//If we want to approximate the integral value for a function, we replace "ptcount++" with the function evaluation.
		}
	}
	return 0;
}


int main(int argc,char *argv[]){

int n=atoi(argv[1]);
int d=(1<<n);
//n,d are set by using argv[1].

int logN=atoi(argv[2]);
//logN denotes log2(N) set by argv[2].

int i;
for( i=1;i<=d;i++){
	int index=0;
	int parity=i;
	int flag=0;
	while( flag!=1){
		if(0==parity%2){
			parity=parity/2; index++;
		}else{flag=1;}
	}
	decomp[i][0]=index;
	decomp[i][1]=parity;
}
//Filling the array of decomp[][]

int index;
for( index=1;index<=n+1;index++){
	int j;
	for( j=0;j<(1<<(index-1));j++){
		double seq=2*cos(M_PI*(2*sigma(index,j+1)-1)/(double)(4*(1<<index-1)));
		Cos[index-1][j+1]=seq;
		Cosinv[index-1][j+1]=1/seq;
	}
}
//Filling the arrays Cos[][] and Cosinv[][]

double disc=pow(2.0*d,d/2.0)/sqrt(2.0);
//Calculating the discriminant of the Vandermonde matrix of the Chebychev roots

double N=pow(2.0,(logN));
double L=pow((N*disc),1.0/(double)d)/2.0; //Calculate the scaling number s(N)

for( i=1;i<=d;i++){
	//Initializing the vectors beta[n] and gamma[n]
	b[n][i]=-L;
	c[n][i]=L;
}

int j;
for( j=n-1;j>=0;j--){
	//Calculating the initial lower and upper bounds defining the point set inside the unit cube
	int h;
	for( h=1;h<=(1<<j);h++){
		b[j][h]=(b[j+1][h]+b[j+1][h+(1<<j)])/2.0;
		c[j][h]=(c[j+1][h]+c[j+1][h+(1<<j)])/2.0;
	}
}

ptcount=0; //Initializing the variable "ptcount"

Enum(1,n); //Counting the number of the Chebyshev-Frolov lattice points

printf("dim=%d, log2(N)=%d, nodes=%lld\n",d,logN,ptcount); //This is the output.

return 0;

}