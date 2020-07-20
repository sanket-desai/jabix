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

import java.util.*;
import java.io.*;

public class GenomicInterval {
  String chromosome="";
  long intervalstart=0, intervalend=0;

  public GenomicInterval(){}
  public GenomicInterval(String interval) throws GenomicIntervalException { //format chrom:start-end
    if ( interval.indexOf(':') == -1 ||  interval.indexOf('-') == -1) {
      throw new GenomicIntervalException();
    }
    this.chromosome=interval.substring(0, interval.indexOf(':'));
    this.intervalstart=Long.parseLong( interval.substring(interval.indexOf(':')+1, interval.indexOf('-')) );
    this.intervalend=Long.parseLong( interval.substring(interval.indexOf('-')+1 ) );
    if (this.intervalstart > this.intervalend ) {
      throw new GenomicIntervalException("End position should be greater than start position in the genomic interval!!");
    }
  }
  public GenomicInterval(String chr, long st, long en) throws GenomicIntervalException {
    this.chromosome=chr;
    this.intervalstart=st;
    this.intervalend=en;
    if (this.intervalstart > this.intervalend ) {
      throw new GenomicIntervalException("End position should be greater than start position in the genomic interval!!");
    }
  }
  long getGenomicIntervalStart(){
    return this.intervalstart;
  }
  long getGenomicIntervalEnd(){
    return this.intervalend;
  }
  String getGenomicIntervalChromosome(){
    return this.chromosome;
  }
  void setGenomicIntervalStart(long gest){
    this.intervalstart=gest;
  }
  void setGenomicIntervalEnd(long geen){
    this.intervalend=geen;
  }
  void setGenomicIntervalChromosome(String gech){
    this.chromosome=gech;
  }
  boolean isWithinRange(String chr, long pos){
    boolean inrange=false;
    if (this.chromosome.equals(chr) && pos >= this.intervalstart && pos <= this.intervalend ) {
        inrange=true;
    }
    return inrange;
  }
  GenomicInterval getOverlap(GenomicInterval gi){
    GenomicInterval geo=new GenomicInterval();
    //Condition 1: if this interval spans within the argument genomic interval
    if (gi.getGenomicIntervalChromosome().equals(this.getGenomicIntervalChromosome()) && gi.getGenomicIntervalStart() <= this.getGenomicIntervalStart() && gi.getGenomicIntervalEnd() >= this.getGenomicIntervalEnd() ) {
      for (long i=gi.getGenomicIntervalStart(); i<=gi.getGenomicIntervalEnd() ; i++ ) {
        if (i == this.getGenomicIntervalStart() ) {
          geo.setGenomicIntervalChromosome(gi.getGenomicIntervalChromosome());
          geo.setGenomicIntervalStart(i);
        }
        if (i == this.getGenomicIntervalEnd()) {
          geo.setGenomicIntervalEnd(i);
          break;
        }
      }
    }
    //Condition 2: if argument interval is spanning within this interval
    else if (gi.getGenomicIntervalChromosome().equals(this.getGenomicIntervalChromosome()) && gi.getGenomicIntervalStart() >= this.getGenomicIntervalStart() && gi.getGenomicIntervalEnd() <= this.getGenomicIntervalEnd() ) {
      for (long i=this.getGenomicIntervalStart(); i <= this.getGenomicIntervalEnd() ; i++ ) {
        if (i == gi.getGenomicIntervalStart() ) {
          geo.setGenomicIntervalStart(i);
          geo.setGenomicIntervalChromosome(gi.getGenomicIntervalChromosome());
        }
        if (i == gi.getGenomicIntervalEnd()) {
          geo.setGenomicIntervalEnd(i);
          break;
        }
      }
    }
    //Condition 3: if any overlap between two A partly spans B or B partly spans A - both conditions are checked below
    else if (this.hasOverlap(gi)) {
      if (gi.getGenomicIntervalStart() >= this.getGenomicIntervalStart() && gi.getGenomicIntervalEnd() <= this.getGenomicIntervalEnd() ) {
        for (long i=this.getGenomicIntervalStart(); i <= this.getGenomicIntervalEnd() ; i++ ) {
          if ( i == gi.getGenomicIntervalStart() ) {
            geo.setGenomicIntervalStart(i);
            geo.setGenomicIntervalChromosome(this.getGenomicIntervalChromosome());
          }
          geo.setGenomicIntervalEnd(gi.getGenomicIntervalEnd());
        }
      }
      else if (gi.getGenomicIntervalStart() <= this.getGenomicIntervalStart() && gi.getGenomicIntervalEnd() >= this.getGenomicIntervalEnd()) {
        for (long i=gi.getGenomicIntervalStart(); i<= gi.getGenomicIntervalEnd() ; i++ ) {
          if ( i == this.getGenomicIntervalStart()) {
            geo.setGenomicIntervalStart(i);
            geo.setGenomicIntervalChromosome(gi.getGenomicIntervalChromosome());
          }
          geo.setGenomicIntervalEnd(this.getGenomicIntervalEnd());
        }
      }
    }
    return geo;
  }
  public String toString(){
    String genointerval="";
    genointerval=this.getGenomicIntervalChromosome()+":"+ Long.toString(this.getGenomicIntervalStart())+"-"+Long.toString(this.getGenomicIntervalEnd());
    return genointerval;
  }
  //Merge, checks if overlap exists if so returns the merged interval else returns "this"
  public GenomicInterval mergeGenomicInterval(GenomicInterval g){
    GenomicInterval og=this;
    if(this.hasOverlap(g)){
      if (this.getGenomicIntervalStart() <= g.getGenomicIntervalStart()) {
        og.setGenomicIntervalStart(this.getGenomicIntervalStart());
      }
      else if (this.getGenomicIntervalStart() > g.getGenomicIntervalStart()) {
        og.setGenomicIntervalStart(g.getGenomicIntervalStart());
      }
      if (this.getGenomicIntervalEnd() >= g.getGenomicIntervalEnd()) {
        og.setGenomicIntervalEnd(this.getGenomicIntervalEnd());
      }
      else if (this.getGenomicIntervalEnd() < g.getGenomicIntervalEnd()) {
        og.setGenomicIntervalEnd(g.getGenomicIntervalEnd());
      }
    }
    return og;
  }
  boolean hasOverlap(GenomicInterval gi){
    boolean overlap=false;
    if (this.isWithinRange(gi.getGenomicIntervalChromosome(), gi.getGenomicIntervalStart()) || this.isWithinRange(gi.getGenomicIntervalChromosome() , gi.getGenomicIntervalStart()) ) {
      overlap=true;
    }
    else if (gi.getGenomicIntervalChromosome().equals(this.getGenomicIntervalChromosome()) && gi.getGenomicIntervalStart() < this.getGenomicIntervalStart() && gi.getGenomicIntervalEnd() > this.getGenomicIntervalEnd() ) {
      overlap=true;
    }
    return overlap;
  }
  //Returns number of variants within genomic interval; argument: TabixIndexedVCF file (bgzipped)
  int getNumberOfVariants(String bgzvcf) throws IOException {
    TabixReader tr=new TabixReader(bgzvcf);
    TabixReader.Iterator iter = tr.query(this.toString());
    int numvariants=0;
    String s;
    while (iter != null && (s = iter.next()) != null){
      numvariants++;
    }
    return numvariants;
  }
  long getLength(){
    return this.getGenomicIntervalEnd()-this.getGenomicIntervalStart()+1;
  }
  long size(){
    return this.getGenomicIntervalEnd()-this.getGenomicIntervalStart()+1;
  }
  List<VCFReader.Variant> getVariantsFromInterval(String bgzvcf) throws IOException, VCFReader.FormatException {
    List<VCFReader.Variant> vars=new ArrayList<VCFReader.Variant>();
    TabixReader tr=new TabixReader(bgzvcf);
    String s;
    TabixReader.Iterator iter=tr.query(this.toString());
    while(iter != null && ( s=iter.next() ) != null ){
      vars.add(new VCFReader.Variant(s));
    }
    return vars;
  }
  public static class GenomicIntervalException extends Exception{
    public GenomicIntervalException(){
      super("Genomic interval provided must be in the format <chrom:start-end>");
    }
    public GenomicIntervalException(String e){
      super(e);
    }
  }
/*
  public static void main(String[] args){
    try{
      GenomicInterval gi=new GenomicInterval("chr22", 24235200, 24236300);
      List<VCFReader.Variant> vars=new ArrayList<VCFReader.Variant>();
      vars=gi.getVariantsFromInterval(args[0]);
      for(VCFReader.Variant v: vars){
        System.out.println(v.toString());
      }
      System.out.println(gi.toString());
      System.out.println(gi.getNumberOfVariants(args[0]));
    }
    catch(GenomicIntervalException ge){
      ge.printStackTrace();
    }
    catch(VCFReader.FormatException fe){
      fe.printStackTrace();
    }
    catch(IOException ie){
      ie.printStackTrace();
    }
  }
  */
}
