import java.util.Random;

public class Photon 
{
	static int width=1000,height=1000;int intensity=170;
	Random rand;
	static double c=10;
	double lambda;
	double pol;
	int[] color;
	double[]loc;
	String bb="ray";//"reflect";//"torus", 
	double dir;// angle from horizontal
	
	public Photon(double l,double[]lo,double d)
	{
		lambda=l;
		loc=lo;
		dir=d;
		rand=new Random();
		pol=rand.nextDouble();
		color=rgb(l,intensity);
	}
	
	public Photon(double[] lo, double d) 
	{
		loc=lo;
		dir=d;
		rand=new Random();
		pol=rand.nextDouble();
		lambda=400*(1+rand.nextDouble());
		color=rgb(lambda,intensity);
	}

	private static int[] rgb(double l, int intensity2) 
	{
		double lb=400,ub=800;
		int c=(int) ((ub-l)/(ub-lb)*6*intensity2)	;int i=0;;
		int		h=c/intensity2;
			
		int[] col=new int[3];
		while(h>1)
		{
			h-=2;
			c-=2*intensity2;
			i++;
		}
		int j=(i+1)%3;
		if(h==0)
		{
			col[i]=intensity2;
			col[j]=c;
		}
		else
		{
			col[j]=intensity2;
			col[i]=2*intensity2-c;
		}
		//System.out.println("lambda="+l+" leads to ("+col[0]+","+col[1]+","+col[2]+")");
		return col;
	}

	public void move(Medium m)
	{
	
		double[]gradient=m.gradient(loc,lambda);
		double normal=Math.atan2(gradient[1], gradient[0]),
				alpha=(dir-normal+3*Math.PI)%(2*Math.PI)-Math.PI,
				n1=m.nindex((int)loc[0],(int)loc[1],lambda),
				norm=Math.sqrt(Math.pow(gradient[0],2)+Math.pow(gradient[1], 2)),
				n2=n1+norm*Math.signum(Math.cos(alpha)),
				n3=n1+norm*Math.signum(Math.cos(alpha))*c/n1,
				beta=Math.sin(alpha)*n1/n3,
				R=1;	
	
		if(Math.abs(beta)<=1) 
		{
			beta=Math.asin(beta);
			if(Math.abs(alpha)>Math.PI/2)beta=(Math.PI-beta)%(2*Math.PI);
			R=pol*Math.pow((n1*Math.cos(alpha)-n3*Math.cos(beta))/(n1*Math.cos(alpha)+n3*Math.cos(beta)),2)
					+(1-pol)*Math.pow((n3*Math.cos(alpha)-n1*Math.cos(beta))/(n3*Math.cos(alpha)+n1*Math.cos(beta)),2);
		}
		else beta=Math.PI/2*Math.signum(alpha);
		//System.out.println("alpha="+alpha+"beta="+beta);
		if(rand.nextDouble()<R)//reflect
		{	//System.out.print("dir="+dir+", "+"normal="+normal+", ");
			//System.out.println("alpha="+alpha+", beta="+beta+"->R="+R+"->reflect!");
			dir=(2*normal-dir+2*Math.PI)%(2*Math.PI)-Math.PI;
			double speed=Math.max(-height, Math.min(c/n1,height));
			loc[0]+=speed*Math.cos(dir);
			loc[1]+=speed*Math.sin(dir);
		}
		else//refract
		{
			dir=(normal+beta+3*Math.PI)%(2*Math.PI)-Math.PI;
			double speed=c/n2;
			loc[0]+=speed*Math.cos(dir);
			loc[1]+=speed*Math.sin(dir);
		}
		if(bb.equals("ray"))
		{if(loc[0]>=width-1||loc[1]>=height-1||loc[0]<0||loc[1]<0) {double r=rand.nextDouble(c);dir=rand.nextDouble(2*Math.PI);loc[0]=r*Math.cos(dir)*LightShow.start[0];loc[1]=r*Math.sin(dir)+LightShow.start[1];}}
		
	//	System.out.println("loc("+loc[0]+","+loc[1]+")");
		else if(bb.equals("torus"))
		{
			if(loc[0]>=width-1)loc[0]-=width-1;
			else if(loc[0]<0)loc[0]+=width-1.000000000001;
			if(loc[1]>=height-1)loc[1]-=height-1;
			else if(loc[1]<0)loc[1]+=height-1.000000000001;
		}
		
		else //mirror
		{
			if(loc[0]>=width-1){loc[0]=2*width-2-loc[0];dir=(Math.PI-dir);}
			else if(loc[0]<0) {loc[0]*=-1;dir=-(Math.PI+dir);}
			if(loc[1]>=height-1) {loc[1]=2*height-2-loc[1];dir*=-1;}
			else if(loc[1]<0) {loc[1]*=-1;dir*=-1;}
		}
		//System.out.println("loc("+loc[0]+","+loc[1]+")");
	}

	public void draw(int[][][] frame) 
	{
		int x0=(int)loc[0],y0=(int)loc[1];
		for(int i=0;i<3;i++)for(int x=x0;x<x0+2;x++)for(int y=y0;y<y0+2;y++)
		{
			frame[x][y][i]=frame[x][y][i]+color[i];
		//	System.out.println(loc[0]+","+loc[1]+"="+frame[x][y][i]);
		}
		
		
	}

	public static Photon[] radial(int i,int j, int k) 
	{
		Photon[]out=new Photon[k];
		for(int n=0;n<k;n++)
			out[n]=new Photon( new double[] {i,j},2*Math.PI/k*n-Math.PI);
		return out;
	}
	public static Photon[] ray(int i,int j,double dir,int w,int k)
	{
		Photon[]out=new Photon[k];
		double l=400, step=1.0*w/k,lstep=400*step, y=j;
		for(int n=0;n<k;n++)
		{	
			out[n]=new Photon(new double[] {i,y},dir);
			y+=step;
			l+=lstep;
			if(l>800)l-=400;
			
		}
		return out;
	}
	
	
}
