import java.awt.Color;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;
import java.text.DecimalFormat;

import javax.imageio.ImageIO;

public class LightShow 
{
	static int it=10000;
	static Photon[]photons;
	static Medium medium;
	static int width,height;
	static int[]start= {0,0};//
	static String name="rd",//"conway",//"escher",//"techlab",//
			type="png";
	static double blur=0.8;
	
	public static void main(String[] args) 
	{
		DecimalFormat df=new DecimalFormat("00000");
		double a=0.0,b=0.15, r=0.001,phibound=Math.atan2(4,-1),phistart=-Math.PI/4,sign=1,phi=-Math.PI/4,step=0.01;
		String b1;
		//double[][][]ce=Medium.ce;
		medium=Medium.newCA(a*2,b*2,0.25);
		System.out.println("medium done");
		//Medium.linear(0, 0, 1, 0);
				//Medium.fromFile("Image00000.png", 1,0.5,2); 
				//Medium.prismring(4,400,1, 0.5,0.05,Photon.c);
				//Medium.radialRange();
				//Medium.Voronoise(96, Photon.c);
				//Medium.randomDrops(Photon.c); 
				//Medium.interference(1.4, 0.8);
				//Medium.circle(0,1,Math.pow(10, 0));
				//Medium.prisma(new double[] {250, 400}, new double[] {750,400},new  double[] {500, 400+500*Math.sqrt(3)/2},a,b,Photon.c);//
		photons=//Photon.radial(start[0],start[1],2000000);
				Photon.ray(start[0], start[1],Math.PI/4, 10, 500000);
		
		int[][][] frame=medium.setUp();
		width=frame.length;
		height=frame[0].length;
	//	for(int k=0;k<50;k++)	medium.MandelIterate();
	//	printmedium(medium);
		for (int i=0;i<it;i++)
		{			//medium=Medium.linear(0.001*Math.cos(i/240.0), 0.001*Math.sin(i/240.0), 1, 0);

			//phi+=sign*step/r;
			//if(phi>phibound) {r+=(step); phi=phibound;sign*=-1;}
			//else if(phi<phistart) {r+=(step); phi=phistart;sign*=-1;}
			//a=1+r*Math.cos(phi);b=r*Math.sin(phi);
			System.out.println(i);
			for(int j=0;j<photons.length;j++)
			{
				photons[j].move(medium);
				photons[j].draw(frame);
			}
			//if(i%24==0) {medium=Medium.fromFile("Image"+df.format(i/24)+".png", 1,0.2,2);}//printmedium(medium);}
		medium.blur(0.972);
		if(i%24==23) {	 medium.updateCA((i/24)%2,a,b);}
		//	printmedium(medium);}
			//medium=Medium.prismring(4,400,1+i/10.0, .5,.05,Photon.c);Medium.prisma(new double[] {250, 400}, new double[] {750,400},new  double[] {500, 400+500*Math.sqrt(3)/2},a,b,Photon.c);//
			if(b<0)b1="";else b1=" ";
			b1+=df.format(b);
			print(frame,i," "+df.format(a),b1);
			blur(frame);
		}

	}
	private static void blur(int[][][] frame) 
	{
		for(int i=0;i<width;i++)
			for(int j=0;j<height;j++)
				for(int k=0;k<3;k++)
					frame[i][j][k]*=blur;
		
	}
	private static void print(int[][][] frame, int t, String a, String b) 
	{
		int[]c=new int[3];
		BufferedImage image=new BufferedImage(width,height, BufferedImage.TYPE_INT_BGR);
		File file=new File(name+t+"."+type);
		for(int i=0;i<width;i++)
			for(int j=0;j<height;j++)
			{
				for(int k=0;k<3;k++)c[k]=256*frame[i][j][k]/(frame[i][j][k]+256);
				image.setRGB(i, j, new Color(c[0],c[1],c[2]).getRGB());
				//System.out.println(frame[i][j][0]+", "+frame[i][j][1]+", "+frame[i][j][2]);
				//int max=Math.max(frame[i][j][0],Math.max(frame[i][j][1],frame[i][j][2]));
				//if(max>8000000)image.setRGB(i, j, Color.white.getRGB());
			//	else if (max>255)image.setRGB(i,j,new Color(frame[i][j][0]*255/max,frame[i][j][1]*255/max,frame[i][j][2]*255/max).getRGB());
			//	else image.setRGB(i,j,new Color(frame[i][j][0],frame[i][j][1],frame[i][j][2]).getRGB());
			}
	//	Write.write(image, "a="+a, 30, 30, 25);
	//	Write.write(image, "b="+b, 30, 85, 25);
		try {
			ImageIO.write(image, type, file);
		} catch (IOException e) {
			System.out.println("IOException: Problems saving file "+name+t);
			e.printStackTrace();
		}
	}
	private static void printmedium(Medium m)
	{
		BufferedImage image=new BufferedImage(width,height, BufferedImage.TYPE_3BYTE_BGR);
		File file=new File(name+"medium"+"."+type);
		for(int i=0;i<width;i++)
			for(int j=0;j<height;j++)
			{
				//System.out.println("green="+(int) (255*m.nIndex[i][j][1]));
				image.setRGB(i,j,new Color((int) Math.min(Math.max(0,255*(m.nIndex[i][j][0]-1)),255),(int) Math.min(255,Math.max(0,255*(m.nIndex[i][j][1]+0.25))),0).getRGB());
			}
	
		try {
			ImageIO.write(image, type, file);
		} catch (IOException e) {
			System.out.println("IOException: Problems saving file "+name+"medium");
			e.printStackTrace();
		}
	}


}
