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

public class SAMAlignmentRecord{
  String[] salignmentrecord_=new String[] {};
  public SAMAlignmentRecord(){}
  public SAMAlignmentRecord(String se) throws SAMAlignmentRecordException {
    this.salignmentrecord_=se.split("\t");
    if (this.salignmentrecord_.length < 11 ) {
      throw new SAMAlignmentRecordException("Number of fields are less in the SAM Alignment Record. Please refer SAM format documentation!!");
    }
  }
  String getQueryName(){
    return this.salignmentrecord_[0];
  }
  String getReadName(){
    return this.salignmentrecord_[0];
  }
  int getFlag(){
    return Integer.parseInt(this.salignmentrecord_[1]);
  }
  String getReferenceSequenceName(){
    return this.salignmentrecord_[2];
  }
  int getAlignedPosition(){
    return Integer.parseInt(this.salignmentrecord_[3]);
  }
  int getMappingQuality(){
    return Integer.parseInt( this.salignmentrecord_[4] );
  }
  String getCIGARString(){
    return this.salignmentrecord_[5];
  }
  String getMateReadName(){
    return this.salignmentrecord_[6];
  }
  int getMateAlignedPosition(){
    return Integer.parseInt( this.salignmentrecord_[7] );
  }
  int getTemplateLength(){
    return Integer.parseInt(this.salignmentrecord_[8]);
  }
  String getReadSequence(){
    return this.salignmentrecord_[9];
  }
  Integer getReadSequenceLength(){
    return this.salignmentrecord_[9].length();
  }
  String getQualityString(){
    return this.salignmentrecord_[10];
  }
  String getOptionalFiledAsString(){
    String optionalfield="";
    if (salignmentrecord_.length==12) {
      optionalfield=this.salignmentrecord_[11];
    }
    return optionalfield;
  }
  public String toString(){
    return String.join("\t",this.salignmentrecord_);
  }
  public class SAMAlignmentRecordException extends Exception {
    public SAMAlignmentRecordException(){ super("SAM Record format exception. Please refer http://samtools.github.io/hts-specs/SAMv1.pdf ."); }
    public SAMAlignmentRecordException(String er){ super(er); }
  }
}
