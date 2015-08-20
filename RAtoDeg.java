/**
   author: Larry Li
   date: June 4th 2015

   This is a small program that converts RA coordinates in hour/min/sec format
   to degrees for the GALFA cubes - since they're all in degrees

*/

import java.util.Scanner;

public class RAtoDeg {
    public static final void main (String[]args){
	//fields
	double inHour;
	double inMin;
	double inSec;
	Scanner in = new Scanner(System.in);
	double outDegFull;

	System.out.println("Input hour");
	inHour = in.nextDouble();
	System.out.println("Input min");
	inMin = in.nextDouble();
	System.out.println("Input sec");
	inSec = in.nextDouble();

	// from hms to degrees
	outDegFull = inHour * 15;
	outDegFull += inMin * 0.25;
	outDegFull += inSec * (15.0/3600);

	System.out.println();
	System.out.println((int)inHour + "h " + (int)inMin + "m " + inSec + "s is");
	System.out.println(outDegFull + " degrees");
	
	double outDeg;
	double outArcMin;
	double outArcSec;

	// from hms to degrees arcMin arcSec

	outDeg = inHour * 15;
	outArcMin = inMin * 15;
	outArcSec = inSec * 15;
	if (outArcSec >= 60){
	    outArcMin += ((int)outArcSec) / 60;
	    outArcSec = outArcSec % 60;
	}
	if (outArcMin >= 60){
	    outDeg += ((int)outArcMin) / 60;
	    outArcMin = outArcMin % 60;
	}

	//degree sign
	final String DEG = "\u00b0";
	
	System.out.println((int)outDeg + DEG + " " + (int)outArcMin + "\' " 
			   + outArcSec + "\"");
	System.out.println("------------------------------------------------");
	
	
	
    }
}
