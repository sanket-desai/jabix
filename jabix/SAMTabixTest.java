package jabix;


public class SAMTabixTest{

  public static void main(String[] args) {
    if (args.length < 1) {
      System.out.println("Usage: java -cp .:sam.jar TabixReader <in.gz> [region]");
      System.exit(1);
    }
    try {
      TabixReader tr = new TabixReader(args[0]);
      String s;
      TabixReader.Iterator iter = tr.query(args[1]); // get the iterator
      while (iter != null && (s = iter.next()) != null){
          //System.out.println(s);
          SAMAlignmentRecord srec=new SAMAlignmentRecord(s);
          System.out.println(srec.getReferenceSequenceName()+" "+srec.getAlignedPosition()+" "+srec.getMateReadName()+" "+ srec.getMateAlignedPosition()+" "+srec.getTemplateLength());
        }
      }
    catch (Exception e) {
      e.printStackTrace();
    }
  }
  //Using fastqToPhred from net.sf.samtools.SAMUtils
}
