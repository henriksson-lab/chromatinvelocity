package chromvel;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.Set;
import java.util.TreeMap;
import java.util.zip.GZIPInputStream;

/**
 * 
 * Generate atac_fragments.tsv from a chromatin velocity BAM input file
 * 
 * Afterwards....
 * compress with: bgzip foo.tsv
 * then index with: tabix -p vcf foo.tsv.gz
 * 
 * Should aim to get 		
 * tnH.atac_fragments.tsv.gz
 * tn5.atac_fragments.tsv.gz
 * 
 * @author Johan Henriksson
 *
 */
public class BamToFragment {
	

	public static BufferedReader openBR(File fInput) throws IOException {
		InputStream is=new FileInputStream(fInput);
		if(fInput.getName().endsWith(".gz"))
			is = new GZIPInputStream(is);
		BufferedReader br = new BufferedReader(new InputStreamReader(is));
		return br;
	}
	
	public static void read(
			File fBAM, 
			TreeMap<String,String> mapReadBarcode,
			File fOut)
					throws IOException {
		
		PrintWriter pw=new PrintWriter(fOut);
		
		TreeMap<String, Integer> lastPosForBC=new TreeMap<>();
		
		Process p = Runtime.getRuntime().exec("samtools view "+fBAM);
		BufferedReader inp = new BufferedReader( new InputStreamReader(p.getInputStream()) );

		//Example line:
		//A00689:445:HNTW5DRXY:1:1114:21359:1125	99	chr1	9996	0	50M	=	10179	232	TCCCATAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAA	F:F,FF,FF:FFFF:FFF:FFFFFFFFFF:FF,FFFFF,FFFFFF,FFFF	NM:i:1	MD:Z:3G46	AS:i:46	XS:i:46	CR:Z:GCAAGGTGTTAACTGC	CY:Z:FFFFFFFFFFFFFFFF	CB:Z:CGGCTCACAACTAGAA-1	RG:Z:lib6:MissingLibrary:1:HNTW5DRXY:1
		
		int dedupCount=0;
		int ignoredCount=0;
		
		String line;
		int readRecords=0;
		while((line=inp.readLine())!=null) {
			//Update user about progress
			readRecords++;
			if(readRecords%1000000 == 0){
				//Calculate progress
				System.out.println("records so far: "+readRecords);
				System.out.flush();
			}
			
			String[] parts=line.split("\t", 0);
			String readname=parts[0];
			String read_chr=parts[2];
			String read_from=parts[3];
			String read_tlen=parts[8];
			
			//System.out.println("*"+readname+"*");
			
			/*
			System.out.println(mapReadBarcode.keySet().iterator().next());
			System.out.println("@"+readname);
			System.out.println();*/
			
			String cb=mapReadBarcode.get("@"+readname);
			if(cb!=null) {

				try {
					Integer lastpos=lastPosForBC.get(cb);
					
					int read_from_int=Integer.parseInt(read_from); //Will this always work?
					if(lastpos==null || lastpos!=read_from_int) {
						
						pw.println(
								read_chr+"\t"+
								read_from+"\t"+
								(Integer.parseInt(read_from)+Integer.parseInt(read_tlen))+"\t"+
								cb+"\t"+
								"1"
						);
						
						lastPosForBC.put(cb, read_from_int);
					} else {
						//System.out.println("Deduplicated read for bc "+cb);
						dedupCount++;
					}
					
				} catch (NumberFormatException e) {
					//System.out.println("Ignored");
					ignoredCount++;
				}
			} else {
				//System.out.println("Missing bc for read: "+ readname);
			}

		}
		inp.close();
		pw.close();
		
		System.out.println("total records read: "+readRecords+" ignored:"+ignoredCount+" dedup:"+dedupCount);
	}
	
	/**
	 * 
	 * 
	 * @param bcfile
	 * @return
	 * @throws IOException
	 */
	public static TreeMap<String,Integer> countFragmentPerBC(File bcfile) throws IOException{
		
		TreeMap<String,Integer> themap=new TreeMap<>();
		
		BufferedReader br=openBR(bcfile);
		
		String line;
		while((line=br.readLine())!=null) {
			line.substring(0,line.indexOf(" "));
			String bc=br.readLine();
			br.readLine();
			br.readLine();
			
			Integer prevcnt=themap.get(bc);
			if(prevcnt==null)
				themap.put(bc, 1);
			else
				themap.put(bc, 1+prevcnt);
		}
		br.close();
		
		return themap;
	}
	
	
	/**
	 * 
	 * 
	 * 
	 * @param bcfile
	 * @param keepbc
	 * @return
	 * @throws IOException
	 */
	public static TreeMap<String,String> readReadBarcodeMap(File bcfile, Set<String> keepbc) throws IOException{

		TreeMap<String,String> themap=new TreeMap<>();		
		BufferedReader br=openBR(bcfile);
		
		String line;
		while((line=br.readLine())!=null) {
			String readname=line.substring(0,line.indexOf(" "));
			String bc=br.readLine();
			br.readLine();
			br.readLine();
			
			if(keepbc.contains(bc)) {
				bc=bc.intern();  //Force string into string pool so it is reused by all barcodes; should save tons of memory
				themap.put(readname, bc);
			}
		}
		br.close();
		
		return themap;
	}
	
	/**
	 * 
	 * 
	 * @param bccount
	 * @param numcells
	 * @return
	 */
	public static TreeMap<String,Integer> filterMinCount(TreeMap<String,Integer> bccount, int numcells) {
		//Count occurrences
		ArrayList<Integer> countlist=new ArrayList<>(bccount.size());
		for(Integer i:bccount.values())
			countlist.add(i);

		//Find cutoff level
		Collections.sort(countlist);
		int cutoff;
		if(countlist.size()<numcells) {
			cutoff=countlist.get(0);
		} else {
			cutoff=countlist.get(countlist.size()-numcells);
		}

		int maxcount=countlist.get(countlist.size()-1);
		System.out.println("Max count for a bc: "+maxcount);
				
		//Filter
		TreeMap<String,Integer> keepbc=new TreeMap<>();
		for(String bc:bccount.keySet()) {
			int cnt=bccount.get(bc);
			if(cnt>=cutoff)
				keepbc.put(bc,cnt);
			if(cnt==maxcount)
				System.out.println(bc);
		}
		System.out.println("------");
		
		return keepbc;
	}
	
	
	
	/**
	 * 
	 * 
	 * 
	 * @param args
	 * @throws IOException
	 */
	public static void main(String[] args) throws IOException {
		
		File readfile=new File("/home/mahogny/temp/chromvel/8_7w_scGETseq_CCTGAGAT.bam");
		File bcfile=new File("/home/mahogny/temp/chromvel/8_7w_scGETseq_R2.fastq.gz");
		File fOut=new File("/home/mahogny/temp/chromvel/8_7w_scGETseq_fragments.tsv");
		
		if(args.length==3) {
			readfile=new File(args[0]);
			bcfile=new File(args[1]);
			fOut=new File(args[2]);
		} else {
			System.out.println("Arguments: aligned.bam  the_R2.fast.gz  out_fragments.tsv");
		}
		
		System.out.println("Count bc");
		TreeMap<String,Integer> bccount = countFragmentPerBC(bcfile);
		//System.out.println(bccount);
		
		System.out.println("Filter bc");
		bccount = filterMinCount(bccount, 8000);
		//System.out.println(bccount);

		System.out.println("Construct read-to-bc map");
		TreeMap<String,String> readbcmap = readReadBarcodeMap(bcfile, bccount.keySet());
		

		System.out.println("Number of reads in map: "+readbcmap.size());
		
		System.out.println("Generate fragment file");
		read(readfile, readbcmap, fOut);
	
		System.out.println("Done");
	}

}
