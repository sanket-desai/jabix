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

public class GTFRecordSAMFeature{
  private List<SAMAlignmentRecord> featuresamrecords;
  private GTFReader.GTFRecord grec;
  private String tabixsam;
  public GTFRecordSAMFeature(){}
  //vcftab -> tabix indexed vcf . Input should be a VCF file
  public GTFRecordSAMFeature(GTFReader.GTFRecord rec, String samtab) throws IOException, VCFReader.FormatException, SAMAlignmentRecord.SAMAlignmentRecordException {
    TabixReader tr=new TabixReader(samtab+".gz");
    String qu=rec.getChrom()+":"+Integer.toString(rec.getStart())+"-"+Integer.toString(rec.getEnd());
    TabixReader.Iterator iter=tr.query(qu);
    String s;
    while(tr != null && ( s=iter.next()) != null ){
      this.featuresamrecords.add(new SAMAlignmentRecord(s));
    }
    this.grec=rec;
    this.tabixsam=samtab;
  }
  SAMAlignmentRecord getGTFRecordSAMAlignmentRecord(int i) throws IndexOutOfBoundsException {
    return featuresamrecords.get(i);
  }
  List<SAMAlignmentRecord> getGTFRecordVariants(){
    return featuresamrecords;
  }
  int getNumberOfSAMAlignmentRecords(){
    return featuresamrecords.size();
  }

}
