import java.awt.Color;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;
import java.util.Random;

import javax.imageio.ImageIO;

public class Medium 
{
	static int width=1000, height=1000;
	double[][][] nIndex;
	//static double[][][]ce=new double[width][height][2],z=new double[width][height][2];//n=a+b/lambda^2
	static int[][][]ca=new int[100][100][2];
	public Medium(double[][][] n)
	{
		nIndex=n;
	}
	public Medium()
	{
		nIndex=new double[width][height][2];
		for(int i=0;i<width;i++)
			for(int j=0;j<height;j++)
				for(int k=0;k<2;k++)
					nIndex[i][j][k]=1-k;
	}

	public double[] gradient(double[] loc,double lambda) 
	{
		int x=(int)loc[0],
				y=(int)loc[1];
		double n=nindex(x,y,lambda);
		double[]out=new double[] {nindex(x+1,y,lambda)-n,nindex(x,y+1,lambda)-n};
		
		return out;
	}

	double nindex(int x, int y, double lambda) 
	{
		
		return nIndex[x][y][0]+nIndex[x][y][1]/lambda/lambda*160000;
	}

	public int[][][] setUp() 
	{
		int[][][]out=new int[width][height][3];
		
		return out;
	}
	public static Medium interference(double a, double b)
	{
		double[][][]n=new double[2560][1440][2];
		for(int i=0;i<2560;i++)
			for(int j=0;j<1440;j++)
			{
				double d0=Math.sqrt(Math.pow(720-i, 2)+Math.pow(720-j, 2)),
						d1=Math.sqrt(Math.pow(1840-i, 2)+Math.pow(720-j, 2));
			//	if(d<300)
				{
					//double value=(1+Math.sin(d0/20))/2;//Math.sqrt(250000-d*d)/100-4;
					n[i][j][0]=1+(a-1)*((1+Math.sin(d0/25))/2)/( Math.log(d0/20+Math.E));//1+(a-1)*value;//
					n[i][j][1]=b*(1+Math.sin(d1/12))/2/( Math.log(d1/15+Math.E));//(0.5+0.5*Math.sin(0.5+d/45));
				}
				/*else
				{
					n[i][j][0]=1;
					n[i][j][1]=0;
				}*/

			}
		return new Medium(n);
	}

	public static Medium circle(double t,double a, double b) 
	{
		double[][][]n=new double[1080][1080][2];
		for(int i=0;i<1080;i++)
			for(int j=0;j<1080;j++)
			{
				double d=Math.sqrt(Math.pow(540-i, 2)+Math.pow(540-j, 2));
			//	if(d<300)
				{
					double value=(1+Math.sin(t/100+d/25))/2;//Math.sqrt(250000-d*d)/100-4;
					n[i][j][0]=1+(a-1)*(value);//1+(a-1)*value;//
					n[i][j][1]=b*value;//(0.5+0.5*Math.sin(0.5+d/45));
				}
				/*else
				{
					n[i][j][0]=1;
					n[i][j][1]=0;
				}*/

			}
		return new Medium(n);
	}
	public static Medium newCA(double a, double b, double d)
	{
		Random rand=new Random();
		Medium m=new Medium();
		double r=width/2.0/ca.length;
		double[]center=new double[2];
		System.out.println("set up ca");
		for(int i=0;i<ca.length;i++)
		{
			center[0]=r*(1+2*i);
			for(int j=0;j<ca[0].length;j++)
				if(rand.nextDouble()<d) 
				{
					ca[i][j][0]=1;
					center[1]=r*(1+2*j);
					m.addDisk(center,r,a,b,Photon.c);
				}
		}
		return m;
	}
	public void addDisk(double[]center,double r, double a, double b, double c)
	{
		double d,min;
		for(int i=0;i<width;i++)
			for(int j=0;j<height;j++)
			{
				d=Math.sqrt(Math.pow(center[0]-i, 2)+Math.pow(center[1]-j, 2));
				if(d<r)
				{
					min=Math.min(r-d, c);
					nIndex[i][j][0]+=a/c*min;
					nIndex[i][j][1]+=b/c*min;
				}
			}
	}
	public static Medium randomDrops(double c)
	{
		Medium m=new Medium();
		Random rand=new Random();
		int n=rand.nextInt(100);
		double r,x,y,a,b;
		for(int i=0;i<n;i++)
		{
			r=rand.nextDouble(height/2)*rand.nextDouble();
			x=rand.nextDouble(r, width-r);
			y=rand.nextDouble(r, height-r);
			a=rand.nextDouble()*rand.nextDouble()*rand.nextDouble();
			b=rand.nextDouble()*rand.nextDouble()*rand.nextDouble()*rand.nextDouble();
			if(a<rand.nextDouble()&&rand.nextBoolean())a*=-1;
			if(b<rand.nextDouble()&&rand.nextBoolean())b*=-1;
			m.addDisk(new double[] {x,y}, r, a, b, c);
		}
		return m;
	}
	public void addPrism(double[]A,double[]B,double[]C,double a,double b,double c)
	{
		double[]AB= vec(A,B),
				BC= vec(B,C),
				CA= vec(C,A);
		double ab=cross(AB,BC),
				bc=cross(BC,CA),
				ca=cross(CA,AB);
		for(int i=0; i<width;i++) {for(int j=0;j<height;j++)
		{
			double distc=cross(AB,new double[] {i-A[0],j-A[1]}),
					dista=cross(BC,new double[] {i-B[0],j-B[1]}),
					distb=cross(CA,new double[] {i-C[0],j-C[1]});
			if(Math.signum(ab)==Math.signum(distc)
					&&Math.signum(bc)==Math.signum(dista)
					&&Math.signum(ca)==Math.signum(distb))
					{
						double value=Math.min(c, Math.min(Math.abs(dista),Math.min(Math.abs(distb), Math.abs(distc))))/c;
						nIndex[i][j][0]+=value*a;
						nIndex[i][j][1]+=value*b;
					}
			
			//System.out.print(n[i][j]);
			}//System.out.println();
		}
	}
	public static Medium prismring(int n, double r1, double r2, double a, double b, double c)
	{
		Medium med=new Medium();
		double[]A=new double[2], B=new double[2], C=new double[2];
		double angle=2*Math.PI/n,phi=angle/2,psi=0;
		if(n%2!=0)phi=0;psi=-angle/2;
		for(int k=0;k<n;k++)
		{
			A[0]=r1*Math.cos(phi)+width/2;
			A[1]=r1*Math.sin(phi)+height/2;
			phi+=angle;
			B[0]=r2*Math.cos(psi)+width/2;
			B[1]=r2*Math.sin(psi)+height/2;
			psi+=angle;
			C[0]=r2*Math.cos(psi)+width/2;
			C[1]=r2*Math.sin(psi)+height/2;
			
			med.addPrism(A, B, C, a, b, c);
		}
		return med;
	}
	public static Medium fromFile(String file, int ir,double rmax, int iphi)
	{
		File source=new File(file);
		double[][][]nin=new double[width][height][2];
		try 
		{
			BufferedImage image=ImageIO.read(source);
			Color c;int[]col;
			int k=(width/2-image.getWidth()/2),l=(height/2-image.getHeight()/2),i1,j1;
			System.out.println("k="+k+", l="+l);
			double r,phi,phical=(Math.atan2(4, -1)+Math.PI/4)/255,a,b;
			for(int i=0;i<width;i++)for(int j=0;j<height;j++)
			{
				i1=Math.max(Math.min(i, width-k-1), k);j1=Math.max(Math.min(j, height-l-1), l);
			//	System.out.println("i="+i1+", j="+j1);
				c=new Color(image.getRGB(i1-k, j1-l));
				col=new int[] {c.getRed(),c.getGreen(), c.getBlue()};
				//r=rmax/255*col[ir];
			//	phi=col[iphi]*phical-Math.PI/4;
				a=rmax*col[ir]/255;
				b=rmax*col[iphi]/255;
				nin[i][j][0]=1+a;//r*Math.cos(phi);
				nin[i][j][1]=b;//r*Math.sin(phi);
			//	System.out.println(nin[i][j][0]+","+nin[i][j][1]);
			}
			
		}
		catch(IOException e){}
		return new Medium(nin);
	}
	public static Medium linear(double mx, double my, double qx, double qy)
	{
		double[][][]nin=new double[width][height][2];
		int midx=height/8, midy=height/2;
		double a,b,move;
		for(int i=0;i<width;i++)
			for(int j=0;j<height;j++)
			{
				a=mx*(i-midx)-my*(j-midy)+qx;
				b=mx*(j-midy)+my*(i-midx)+qy;
				if(b<4-4*a)
				{
					move=(4-4*a-b)/17.0;
					a=a+4*move;
					b=b+move;
				}
				if(b<1-a)
				{
					move=(1-a-b)/2.0;
					a=a+4*move;
					b=b+move;
				}
				nin[i][j][0]=a;
				nin[i][j][1]=b;
			}
		return new Medium(nin);
	}
	public void blur(double b)
	{
		int d=width/ca.length;
		for(int i=d;i<width-d;i++)
			for(int j=d;j<height-d;j++)
			{
				nIndex[i][j][0]*=b;
				nIndex[i][j][0]+=(1-b);
				nIndex[i][j][1]*=b;
			}
	}
	public static Medium prisma(double[]A, double[]B, double[]C, double a, double b, double c)//c=speed of light, diff= difference in refr index
	{
		//double b=64000;//color refr strength
		double[][][]n=new double[width][height][2];
		double[]AB= vec(A,B),
				BC= vec(B,C),
				CA= vec(C,A);
		double ab=cross(AB,BC),
				bc=cross(BC,CA),
				ca=cross(CA,AB);
		for(int i=0; i<width;i++) {for(int j=0;j<height;j++)
		{
			double distc=cross(AB,new double[] {i-A[0],j-A[1]}),
					dista=cross(BC,new double[] {i-B[0],j-B[1]}),
					distb=cross(CA,new double[] {i-C[0],j-C[1]});
			if(Math.signum(ab)==Math.signum(distc)
					&&Math.signum(bc)==Math.signum(dista)
					&&Math.signum(ca)==Math.signum(distb))
					{
						double value=Math.min(c, Math.min(Math.abs(dista),Math.min(Math.abs(distb), Math.abs(distc))))/c;
						n[i][j][0]=1+value*(a-1);
						n[i][j][1]=value*b;
					}
			else { n[i][j][0]=1;n[i][j][1]=0;}
			//System.out.print(n[i][j]);
			}//System.out.println();
		}
		return new Medium(n);
	}
	
	public static double[]vec(double[]A, double[] B)
	{
		double norm=0;
		double[] out= new double[A.length];
		for(int i=0;i<A.length;i++)
		{
			out[i]=B[i]-A[i];
			norm+=out[i]*out[i];
		}
		
		norm=Math.sqrt(norm);
		for(int i=0;i<A.length;i++)
		{
			out[i]/=norm;
		}
		return out;
	}
	
	public static Medium Voronoise(int n,double c)
	{
		Random rand=new Random();
		double[][][]nin=new double[width][height][2];
		double[][] cell=new double[n][5];//x,y,strength, a,b
		double r,phi;
		for(int i=0;i<n;i++)
		{
			cell[i][0]=rand.nextDouble(width);
			cell[i][1]=rand.nextDouble(height);
			cell[i][2]=rand.nextDouble()*0.9+1;	
			r=rand.nextDouble(2)*rand.nextDouble();
			phi=rand.nextDouble(-Math.PI/4,Math.atan2(4, -1));
			
			cell[i][3]=r*Math.cos(phi)+1;
			cell[i][4]=r*Math.sin(phi);		
			System.out.println("center"+cell[i][0]+","+cell[i][1]+", strength"+cell[i][2]+", a"+cell[i][3]+", b"+cell[i][4]);
		}
		int next,second=-1;
		double dist,dist2,d,f;
		for(int x=0;x<width;x++)
			for(int y=0;y<height;y++)
			{
				next=0;dist=dist(x,y,cell[0][0],cell[0][1])*cell[0][2];dist2=width;
				for(int i=1;i<n;i++)
				{
					d=dist(x,y,cell[i][0],cell[i][1])*cell[i][2];
					if(d<dist2)
					{
						if (d<dist)
						{
							second=next;
							next=i;
							dist2=dist;
							dist=d;
						}
						else
						{
							second=i;
							dist2=d;
						}
					}
				}
				f=Math.min(c, (dist2-dist));
				if(f==c)
				{
					nin[x][y][0]=cell[next][3];
					nin[x][y][1]=cell[next][4];
				}
				else
				{
					f=f/c/2+0.5;
					
					nin[x][y][0]=f*cell[next][3]+(1-f)*cell[second][3];
					nin[x][y][1]=f*cell[next][4]+(1-f)*cell[second][4];
				}
				//if(next!=0)System.out.print(next);
			}
		return new Medium(nin);
	}
	public static Medium radialRange()
	{
		double r,phi;
		double[][][]n=new double[width][height][2];
		for(int i=0;i<width;i++)
			for(int j=0;j<height;j++)
			{
				r=dist(i,j,width/2,height/2)*3/height;
				phi=Math.atan2(j-height/2,i-width/2);
				ce[i][j][0]=r*Math.cos(phi);
				ce[i][j][1]=r*Math.sin(phi);
				z[i][j][0]=ce[i][j][0];
				z[i][j][1]=ce[i][j][1];
			//	r=1/r;
				phi=(phi/(2*Math.PI)+0.5)*(Math.PI/4+Math.atan2(4,-1))-Math.PI/4;
				n[i][j][0]=Math.max(0, 1-r/2)*Math.cos(phi)+1;//r*Math.cos(phi)+1;
				n[i][j][1]=Math.max(0, 1-r/2)*Math.sin(phi);//r*Math.sin(phi);
			}
		return new Medium(n);
	}
	public void MandelIterate()
	{
		double r,phi,x,y;
		for(int i=0;i<width;i++)
			for(int j=0;j<height;j++)
				if(!(nIndex[i][j][0]==1&&nIndex[i][j][1]==0))
			{
				
			/*	x=nIndex[i][j][0]-1-ce[i][j][0];
				y=nIndex[i][j][1]-ce[i][j][1];
				r=dist(0,0,x,y);
				phi=Math.atan2(y,x);
				r=Math.min(Math.sqrt(r),2);
				phi=phi/2;*/
				r=dist(0,0,z[i][j][0],z[i][j][1]);
				phi=Math.atan2(z[i][j][1], z[i][j][0]);
				//r=2-2*r;
				//phi=((phi+Math.PI/4)/(Math.PI/4+Math.atan2(4, -1))-0.5)*2*Math.PI;
				r=r*r;
				phi*=2;
				x=r*Math.cos(phi)+(ce[i][j][0]);
				y=r*Math.sin(phi)+ce[i][j][1];
				z[i][j][0]=x;
				z[i][j][1]=y;
				r=dist(0,0,x,y);
				phi=Math.atan2(y, x);
				phi=(phi/(2*Math.PI)+0.5)*(Math.PI/4+Math.atan2(4,-1))-Math.PI/4;
				nIndex[i][j][0]=Math.cos(phi)*Math.max(0, 1-r/2)+1;
				nIndex[i][j][1]=Math.sin(phi)*Math.max(0, 1-r/2);
			}
		
	}
	private static double dist(double x, double y, double d, double e) 
	{
		
		return Math.sqrt(Math.pow(d-x, 2)+Math.pow(e-y, 2));
	}
	public static double cross(double[]v,double[]w)
	{
		return v[0]*w[1]-w[0]*v[1];
	}
	public void updateCA(int par, double a, double b)//par 0 or 1, which is current?
	{
		System.out.println("par="+par);
		double r=width/2.0/ca.length;
		double[]center=new double[2];
		//blur(blur);
		for(int i=1;i<ca.length-1;i++)
		{
			center[0]=r*(1+2*i);
			for(int j=1;j<ca[0].length-1;j++)
			{
				int n=-ca[i][j][par];
				for(int k=-1;k<2;k++)
					for(int l=-1;l<2;l++)
					{
						n+=2*ca[i+k][j+l][par];
					}
				if(n>4&&n<8) 
				{
					ca[i][j][1-par]=1;
					center[1]=r*(1+2*j);
					addDisk(center,r,a,b,Photon.c);
				}
				else ca[i][j][1-par]=0;
			}
		}
		
	}
	
}
