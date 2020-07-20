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

public class GenomicIntervalList {
  List<GenomicInterval> gilist;

  public GenomicIntervalList(){
    gilist=new ArrayList<GenomicInterval>();
  }
  public GenomicIntervalList(List<GenomicInterval> li){
    this.gilist=li;
  }
  public GenomicIntervalList(GTFReader gr) throws GenomicInterval.GenomicIntervalException, IOException, GTFReader.FormatException, GTFParseException {
    this.gilist=new ArrayList<GenomicInterval>();
    while(true){
      try{
        GTFReader.GTFRecord grec=gr.getNextGTFRecord();
        gilist.add(grec.asGenomicInterval());
      }
      catch(NullPointerException ne){
        break;
      }
    }
  }
  public GenomicIntervalList(String fn) throws IOException, GenomicInterval.GenomicIntervalException, GTFReader.FormatException, GTFParseException {
    this.gilist=new ArrayList<GenomicInterval>();
    GTFReader gr=new GTFReader(fn);
    GTFReader.GTFRecord grec=gr.getNextGTFRecord();
    while(true){
      gilist.add(grec.asGenomicInterval());
    }
  }
  GenomicIntervalList getGenomicIntervalList(){
    return this;
  }
  GenomicIntervalList getNonOverlappingGenomicIntervalList(String sizesfi) throws IOException, GenomicInterval.GenomicIntervalException, FileNotFoundException{
    List<GenomicInterval> novgi=new ArrayList<GenomicInterval>();
    BufferedReader sbr=new BufferedReader(new FileReader(sizesfi));
    String s="";
    Map<String, boolean[]> chromboolmap=new HashMap<String, boolean[]>();
    while( sbr != null && (s=sbr.readLine()) != null){
      String[] sl=s.trim().split("\t");
      if(sl.length == 2){
        int csize=Integer.parseInt(sl[1]);
        boolean[] arrabool=new boolean[csize];
        for(int i=0; i < csize; i++){
          arrabool[i]=false;
        }
        chromboolmap.put(sl[0], arrabool);
      }
    }
    for (int ii=0; ii < this.getNumberOfGenomicIntervals() ; ii++ ) {
      GenomicInterval gi=this.gilist.get(ii);
      if (chromboolmap.containsKey(gi.getGenomicIntervalChromosome())) {
        boolean[] boolarr=chromboolmap.get(gi.getGenomicIntervalChromosome());
        for (int win= (int)gi.getGenomicIntervalStart()-1 ; win < gi.getGenomicIntervalEnd() ; win++ ) {
            boolarr[win]=true;
        }
        boolean replacecheck=chromboolmap.replace(gi.getGenomicIntervalChromosome(), chromboolmap.get(gi.getGenomicIntervalChromosome()), boolarr);
        if (!replacecheck) {
          System.out.println("Cannot find Chromosome in the chrom.sizes file! Please use correct file!!");
        }
      }
    }
    Set< Map.Entry< String,boolean[]> > setview = chromboolmap.entrySet();
    for (Map.Entry< String,boolean[]> me:setview){
      int gistart=0, giend=0;
      boolean[] bme= me.getValue();
      boolean featureon=false;
      for(int ind=0; ind < bme.length; ind++ )
      {
        if (bme[ind] && !featureon) {
          gistart=ind+1;
          featureon=true;
        }
        else if ( !bme[ind] && featureon ) {
          giend=ind; //could be simply g++
          featureon=false;
          novgi.add(new GenomicInterval(me.getKey(), gistart, giend));
          gistart=0; giend=0;
        }
        else if (bme[ind] && featureon) {
          giend++;
        }
        else if ( !bme[ind] && !featureon) {}
      }
    }
    return new GenomicIntervalList(novgi);
  }
  int getNumberOfGenomicIntervals(){
    return this.gilist.size();
  }
  int getNumberOfNonOverlappingGenimicIntervals(String chrsifi) throws IOException, GenomicInterval.GenomicIntervalException, FileNotFoundException{
    return this.getNonOverlappingGenomicIntervalList(chrsifi).getNumberOfGenomicIntervals();
  }
  long size(){
    long si=0;
    //GenomicIntervalList gil=this.getNonOverlappingGenomicIntervalList(sizesfi);
    for (int i=0; i < this.getNumberOfGenomicIntervals() ; i++ ) {
      si+=this.getGenomicInterval(i).size();
    }
    return si;
  }
  GenomicInterval getGenomicInterval(int i){
    return this.gilist.get(i);
  }
  void addGenomicInterval(GenomicInterval g){
    this.gilist.add(g);
  }
}
