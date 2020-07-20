/*
   The MIT License

   Copyright (c) 2020 Sanket Desai.

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
 * Email    : desai.sanket12@gmail.com
*/

package jabix;

import java.io.*;
import java.util.*;

public class SAMCoverageMap{
  Map< String, int[] > ref_coveragearray;
  int number_of_alignmentrecords;
  String filename;
  public SAMCoverageMap(){}
  public SAMCoverageMap(String samfi) throws IOException, FileNotFoundException, SAMAlignmentRecord.SAMAlignmentRecordException, SAMHeaderRecord.SAMHeaderRecordException {
    this.ref_coveragearray=new HashMap< String, int[] >();
    this.filename=samfi;
    this.number_of_alignmentrecords=0;
    BufferedReader br=new BufferedReader(new FileReader(this.filename));
    String line="";
    while((line=br.readLine()) != null){
      if(line.startsWith("@")){
        SAMHeaderRecord shr=new SAMHeaderRecord(line);
        if (shr.isHeaderType("SQ")) {
          int seqlength = Integer.parseInt(shr.getHeaderFeatureValue("LN"));
          String seqname= shr.getHeaderFeatureValue("SN");
          //Fill sequname and length into the map.
          ref_coveragearray.put(seqname, new int[seqlength]);
        }
      }
      else{
        SAMAlignmentRecord arec=new SAMAlignmentRecord(line);
        //Module assumes the alignment to be complete / without taking into consideration the the MAPQ / future version should incorporate that!
        try{
          //Integer temparray=new Integer[];
          int[] temparray = this.ref_coveragearray.get(arec.getReferenceSequenceName());
          for (int i = arec.getAlignedPosition()-1 ; i < arec.getAlignedPosition()-1+arec.getReadSequenceLength(); i++ ) {
            temparray[i]=temparray[i]+1;
            this.ref_coveragearray.replace(arec.getReferenceSequenceName(), temparray);
//            this.ref_coveragearray.get(arec.getReferenceSequenceName())[i]=this.ref_coveragearray.get(arec.getReferenceSequenceName())[i]+1;
          }
          //this.ref_coveragearray[arec.getReferenceSequenceName()]=temparray;
        }
        catch(ArrayIndexOutOfBoundsException ie){
          System.out.println("Alignment position out of bounds for the contig reference; alignment position: " +Integer.toString(arec.getAlignedPosition()));
          continue;
        }
        this.number_of_alignmentrecords++;
      }
    }
    ArrayList<String> remref=new ArrayList<String>();
    for (String ke : this.ref_coveragearray.keySet()) {
      int rcsum=0;
      for (int rc : this.ref_coveragearray.get(ke)) {
        rcsum+=rc;
      }
      if (rcsum==0) {
        remref.add(ke);
      }
    }
    for(int ind=0; ind < remref.size(); ind++){
      this.ref_coveragearray.remove(remref.get(ind));
    }
  }
  int getNumberOfAlignmentRecords(){
    return this.number_of_alignmentrecords;
  }
  String getSAMFileName(){
    return this.filename;
  }
  int getNumberOfReferenceContigs(){
    return this.ref_coveragearray.size();
  }
  boolean hasReferenceContig(String refcon){
    return this.ref_coveragearray.containsKey(refcon);
  }
  int getReferenceContigLength(String refcon){
    return this.ref_coveragearray.get(refcon).length;
  }
  String[] getReferenceContigListAsArray(){
    return this.ref_coveragearray.keySet().toArray(new String[this.ref_coveragearray.size()]);
  }
  int[] getCoverage(String refname){
    return this.ref_coveragearray.get(refname);
  }
  double getMeanCoverage(String refname){
    int totalcov=0;
    for (int i =0; i < this.ref_coveragearray.get(refname).length ; i++ ) {
      totalcov=totalcov+this.ref_coveragearray.get(refname)[i];
    }
    return totalcov / this.ref_coveragearray.get(refname).length ;
  }
  double[] getPerKilobaseMeanCoverage(String refname){
    int[] tempcov=this.getCoverage(refname);
    int kblength= (int) Math.ceil((double) tempcov.length / 1000);
    double[] perkbmeancov=new double[kblength];
    int start=0;
    int end=1000;
    for (int i=0; i < kblength; i++) {
      double pkbmean=0;
      int pkbtotrc=0;
      if (end >= tempcov.length) {
        end=tempcov.length;
      }
      for(int rc : Arrays.copyOfRange(tempcov, start, end )){
        pkbtotrc+=rc;
      }
      perkbmeancov[i]= (double) pkbtotrc/(end-start);
      start=start+1000;
      end=end+1000;
    }
    return perkbmeancov;
  }
  int getTotalAlignedReads(String refname){
    int totrc=0;
    for(int ip : this.ref_coveragearray.get(refname)){
      totrc+=ip;
    }
    return totrc;
  }
  //Arguments: comma separated (no space) sample names, comma separated sam files, pathogen name (as per 1070 reference nomenclature), output file
  //In the modified version - Arguments: file name which contains a id separated by complete file path (space separator), pathogen name (as per 1070 reference nomenclature), output file
  public static void main(String[] args){
    try{
      //System.out.println(args[0]);
      Map<String, double[]> sample_pathocov_map=new HashMap<String, double[]>();
      System.out.println(args[0]);
      FileReader fin=new FileReader(args[0]);
      //String[] samplenames=args[0].split(",");
      //String[] filepaths=args[1].split(",");
      BufferedReader br=new BufferedReader(fin);
      List<String> samplenames=new ArrayList<String>();
      List<String> filepaths=new ArrayList<String>();
      String pathoname=args[1];
      String thisLine="";
      try{
      	while ((thisLine = br.readLine()) != null) {
		String[] ss=thisLine.split(" ");
		samplenames.add(ss[0]);
		filepaths.add(ss[1]);
      	}
	}
      catch(Exception e){ e.printStackTrace(); }
      if (samplenames.size() != filepaths.size()) {
        System.out.println("Sample names and sam files provided are not consistent! Please check..");
        System.exit(0);
      }
      for (int filecounter=0; filecounter < samplenames.size() ; filecounter++ ) {
        String samplename=samplenames.get(filecounter);
        SAMCoverageMap scm=new SAMCoverageMap(filepaths.get(filecounter));
        //int [] acov=new int[scm.getReferenceContigLength("Fusobacterium_nucleatum")];
        if (scm.hasReferenceContig(pathoname)) {
          sample_pathocov_map.put( samplename ,scm.getPerKilobaseMeanCoverage(pathoname));
        }
        else{
          System.out.println("Reference contig does not exist!! "+scm.getSAMFileName());
          System.exit(0);
        }
      }
      FileWriter fw=new FileWriter(args[2]);
      for(Map.Entry<String, double[]> me : sample_pathocov_map.entrySet()){
        String mecov=Arrays.toString(me.getValue()).replaceAll("\\s","");
        mecov=mecov.replace("[","");
        mecov=mecov.replace("]","");
        fw.write(me.getKey()+"\t"+ mecov +"\n");
      }
      fw.close();
    }
    catch(IOException ie){
      ie.printStackTrace();
    }
    catch(SAMAlignmentRecord.SAMAlignmentRecordException sre){
      sre.printStackTrace();
    }
    catch(SAMHeaderRecord.SAMHeaderRecordException shre){
      shre.printStackTrace();
    }
    catch(ArrayIndexOutOfBoundsException ae){
      ae.printStackTrace();
    }
  }
}
