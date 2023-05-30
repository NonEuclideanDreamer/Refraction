import java.awt.Color;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;

import javax.imageio.ImageIO;

public class Write 
{
	
	public static int black=Color.black.getRGB(),
			white=Color.white.getRGB();
	static void write(BufferedImage canvas, String text, int x, int y, int size)
	{
		int locx=x,locy=y;
		for(int i=0;i<text.length();i++)
		{
			locx+=writechar(canvas, text.charAt(i),locx,locy,size,white);
		}
	}

	private static int writechar(BufferedImage canvas, char ch, int locx, int locy, int size,int color) 
	{
		//digital display
		/*boolean[] line=new boolean[7];
		if(ch=='0')line=new boolean[]{true,true,true,true,true,true,false};
		else if(ch=='1')line=new boolean[]{false,false,false,true,true,false,false};*/
		if(ch==' ')return size;
		File letterFile=new File(ch+".bmp");
		BufferedImage letter;
		try {
			letter = ImageIO.read(letterFile);	
			double factor=size/400.0;
		for(int i=0;i<letter.getWidth();i++)
			for(int j=0;j<letter.getHeight();j++)
				if(letter.getRGB(i,j)==black)canvas.setRGB((int)(i*factor)+locx, (int)(j*factor)+locy, color);
			return (int) (letter.getWidth()*factor);} 
		catch (IOException e) {
			System.out.println("Can't read "+ch);
			e.printStackTrace();
		}
		
		return 0;
	}
	
}
