/*
   The MIT License

   Copyright (c) 2018 Sanket Desai.

   Permission is hereby granted, free of charge, to any person obtaining
   a copy of this software and associated documentation files (the
   "Software"), to deal in the Software without restriction, including
   without limitation the rights to use, copy, modify, merge, publish,
   distribute, sublicense, and/or sell copies of the Software, and to
   permit persons to whom the Software is furnished to do so, subject to
   the following conditions:

   The above copyright notice and this permission notice shall be
   included in all copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
   NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
   BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
   ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
   CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
   SOFTWARE.
*/

/*
 * Developer: Sanket Desai
 * Email    : desai.sanket12@outlook.com
*/

package jabix;

import java.io.*;
import java.math.*;
import java.util.zip.GZIPInputStream;

public class Utilities{
  //https://stackoverflow.com/a/30507742 - isGZipped
  public static boolean isGZipped(File f) {
    int magic = 0;
    try {
      RandomAccessFile raf = new RandomAccessFile(f, "r");
      magic = raf.read() & 0xff | ((raf.read() << 8) & 0xff00);
      raf.close();
    }
    catch (Throwable e) {
      e.printStackTrace(System.err);
    }
    return magic == GZIPInputStream.GZIP_MAGIC;
  }

  public static boolean isGZipped(String fn) {
    File f=new File(fn);
    int magic = 0;
    try {
      RandomAccessFile raf = new RandomAccessFile(f, "r");
      magic = raf.read() & 0xff | ((raf.read() << 8) & 0xff00);
      raf.close();
    }
    catch (Throwable e) {
      e.printStackTrace(System.err);
    }
    return magic == GZIPInputStream.GZIP_MAGIC;
  }

  public static BigInteger factorial(int j){
    BigInteger fact=BigInteger.valueOf(1);
    for(int i=1;i<=j;i++){
      fact=fact.multiply(BigInteger.valueOf(i));
    }
    return fact;
    /*
    System.out.println("What is I: ");
    System.out.println(i);
    if(i == 0) return 1;      //return 1 when base case is reached
    else return i*fact(i-1);
    */
  }
  //Given mean and count calculate Poisson probability using e~2.7182...
  //mean is expected number of events, count is the observed events
  public static BigDecimal calculatePoissonProbability(double mean, int count){
    //double pdf=0;
    int i=count;
    double e = 2.718281828;
    BigInteger denominator=factorial(i);
    BigDecimal numerator=new BigDecimal(Math.pow(e, -mean) * Math.pow(mean, count));
    //System.out.println(denominator);
    //System.out.println(numerator);
    //double numerator=Math.pow(e, -mean) * Math.pow(mean, count);
    BigDecimal pdf= numerator.divide(new BigDecimal(denominator),5,RoundingMode.HALF_UP);
    return pdf;
  }
  public static double expectedEvents(long sizefeat, long sumfeatsize, int totevents){
    double expevents= ((double) sizefeat / (double) sumfeatsize) * (double) totevents;
    return expevents;
  }
}
