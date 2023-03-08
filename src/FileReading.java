import java.io.*;
import java.math.*;
import java.util.*;

public class FileReading {
		static ArrayList<String> nodes=new ArrayList<>();
		static ArrayList<String> pipes=new ArrayList<>();
		static ArrayList<String> arcs=new ArrayList<>();
		static ArrayList<String> arcsLength=new ArrayList<>();
		static ArrayList<String> elevation=new ArrayList<>();
		static ArrayList<String> pressure=new ArrayList<>();
		static ArrayList<String> demand=new ArrayList<>();
		static ArrayList<String> diameter=new ArrayList<>();
		static ArrayList<String> pipeCost=new ArrayList<>();
		static ArrayList<String> pipeRoughness=new ArrayList<>();
		static String source;
		public static void main(String[] args) throws IOException {
			
			String filePath="F:\\IITB\\MTP\\MTP stage 1\\JalTantra-Code-and-Scripts\\Files\\Data\\Ampl\\d1_Sample_input_cycle_twoloop.dat";
			File file = new File(filePath); 
			Scanner sc=new Scanner(file);
			
			int i=1;
			while(sc.hasNext()){
				if(i==1){
					readnodes(sc);
				}
				else if(i==2){
					readPipes(sc);
				}
				else if(i==3){
					readArcs(sc);
				}
				else if(i==4){
					readDiameter(sc);
				}
				else if(i==5){
					readElevation(sc);
				}
				else if(i==6){
					readPressure(sc);
				}
				else if(i==7){
					readDemand(sc);
				}
				else if(i==8){
					readCost(sc);
				}
				else if(i==9){
					readRoughness(sc);
				}
				else {
					readHead(sc);
					break;
				}
				i++;
			}
			printSetsAndParametersAndModels();
		}
		public static void readnodes(Scanner sc){
			String inp=sc.nextLine();
		    String str[]=inp.split("\\s+");
		    String node="";
		    for(int i=3;i<str.length-2;i++) 
				node=node+str[i]+", ";
			node=node+str[str.length-2];
		    nodes.add(node);
//		    System.out.println(nodes);
		}
		public static void readPipes(Scanner sc){
			String inp=sc.nextLine();
		    String str[]=inp.split("\\s+");
		    String pipe="";
		    for(int i=3;i<str.length-2;i++) 
				pipe=pipe+str[i]+", ";
			pipe=pipe+str[str.length-2];
		    pipes.add(pipe);
//		    System.out.println(pipes);
		}
		public static void readArcs(Scanner sc){
		    String s=sc.nextLine();
			while(sc.hasNext()){
			    s=sc.nextLine();
				if(s.indexOf(";") != -1) 
					break;
			    String str[]=s.split("\\s+");
			    String arc=str[0]+"."+str[1];
			    arcs.add(arc);
			    String len=str[0]+"    ."+str[1]+"    "+str[2];
			    arcsLength.add(len);
			}
			String str[]=s.split("\\s+");
		    String arc=str[0]+"."+str[1];
		    arcs.add(arc);
		    String len=str[0]+"    ."+str[1]+"    "+str[2];
		    arcsLength.add(len);
//		    System.out.println(pipes);
		    sc.nextLine();
		}
		public static void readDiameter(Scanner sc){
			String inp=sc.nextLine();
			while(sc.hasNext()){
			    inp=sc.nextLine();
			    if(inp.indexOf(";") != -1) 
					break;
			    String str[]=inp.split("\\s+");
			    String dia="";
			    for(String s:str) dia=dia+" "+s;
			    diameter.add(dia);
			}
			String str[]=inp.split("\\s+");
		    String dia=str[0]+" "+str[1];
		    diameter.add(dia);
		    sc.nextLine();
		}
		public static void readElevation(Scanner sc){
			String inp=sc.nextLine();
			while(sc.hasNext()){
			    inp=sc.nextLine();
			    if(inp.indexOf(";") != -1) 
					break;
			    String str[]=inp.split("\\s+");
			    String ele="";
			    for(String s:str) ele=ele+" "+s;
			    elevation.add(ele);
			}
			String str[]=inp.split("\\s+");
		    String ele=str[0]+" "+str[1];
		    elevation.add(ele);
		    sc.nextLine();
		}
		public static void readPressure(Scanner sc){
			String inp=sc.nextLine();
			while(sc.hasNext()){
			    inp=sc.nextLine();
			    if(inp.indexOf(";") != -1) 
					break;
			    String str[]=inp.split("\\s+");
			    String pres="";
			    for(String s:str) pres=pres+" "+s;
			    pressure.add(pres);
			}
			String str[]=inp.split("\\s+");
		    String pres=str[0]+" "+str[1];
		    pressure.add(pres);
		    sc.nextLine();
		}
		public static void readDemand(Scanner sc){
			String inp=sc.nextLine();
			while(sc.hasNext()){
			    inp=sc.nextLine();
			    if(inp.indexOf(";") != -1) 
					break;
			    String str[]=inp.split("\\s+");
			    String dem="";
			    for(String s:str) dem=dem+" "+s;
			    demand.add(dem);
			}
			String str[]=inp.split("\\s+");
		    String dem=str[0]+" "+str[1];
		    demand.add(dem);
		    sc.nextLine();
		}
		public static void readCost(Scanner sc){
			String inp=sc.nextLine();
			while(sc.hasNext()){
			    inp=sc.nextLine();
			    if(inp.indexOf(";") != -1) 
					break;
			    String str[]=inp.split("\\s+");
			    String cost="";
			    for(String s:str) cost=cost+" "+s;
			    pipeCost.add(cost);
			}
			String str[]=inp.split("\\s+");
		    String cost=str[0] + " " + str[1];
		    pipeCost.add(cost);
		    sc.nextLine();
		}
		public static void readRoughness(Scanner sc){
			String inp=sc.nextLine();
			while(sc.hasNext()){
			    inp=sc.nextLine();
			    if(inp.indexOf(";") != -1) 
					break;
			    String str[]=inp.split("\\s+");
			    String roughness="";
			    for(String s:str) roughness=roughness+" "+s;
			    pipeRoughness.add(roughness);
			}
			String str[]=inp.split("\\s+");
		    String roughness=str[0]+" "+str[1];
		    pipeRoughness.add(roughness);
		    sc.nextLine();
		}
		public static void readHead(Scanner sc){
			String inp=sc.nextLine();
		    source=inp.split("\\s+")[2];
		}
		public static void printSetsAndParametersAndModels() throws IOException{
			BufferedWriter out = null;
			String gamsFilePath="F:\\IITB\\MTP\\MTP stage 1\\JalTantra-Code-and-Scripts\\Files\\Data\\Gams\\testing.gms";
			
			try {
			    FileWriter fstream = new FileWriter(gamsFilePath, true); //true tells to append data.
			    out = new BufferedWriter(fstream);

				printSets(out);
				printParameters(out);
				printModel(out);
			}

			catch (IOException e) {
			    System.err.println("Error: " + e.getMessage());
			}

			finally {
			    if(out != null) {
			        out.close();
			    }
			}
		}
		static void printSets(BufferedWriter out) throws IOException{
			
			out.write("Sets\n");   //new line
			out.write("\t"+"nodes /");
			for(String str:nodes)
				out.write(str+"");
			out.write("/\n");

			out.write("\t"+"pipes /");
			for(String str:pipes)
				out.write(str+"");
			out.write("/\n");

			out.write("\t"+"src(nodes) /"+source+"/;");
			out.write("\n");
			out.write("alias (src,srcs);\n");
			out.write("alias (nodes,j) ;\n");

			out.write("Set "+"arcs(nodes,j) /");
			for(int i=0;i<arcs.size()-1;i++){
				out.write(arcs.get(i)+", ");
			}
			out.write(arcs.get(arcs.size()-1)+"/\n");
			
			out.write("\n");
		}
		static void printParameters(BufferedWriter out) throws IOException{
			out.write("Parameters\n");
			out.write("\t"+"Len(nodes,j) /");
			for(int i=0;i<arcsLength.size()-1;i++){
				out.write(arcsLength.get(i)+", ");
			}
			out.write(arcsLength.get(arcsLength.size()-1)+"/"+"\n");

			out.write("\t"+"E(nodes) /");
			for(int i=0;i<elevation.size()-1;i++){
				out.write(elevation.get(i)+", ");
			}
			out.write(elevation.get(elevation.size()-1)+"/\n");

			out.write("\t"+"P(nodes) /");
			for(int i=0;i<pressure.size()-1;i++){
				out.write(pressure.get(i)+", ");
			}
			out.write(pressure.get(pressure.size()-1)+"/\n");

			out.write("\t"+"D(nodes) /");
			for(int i=0;i<demand.size()-1;i++){
				out.write(demand.get(i)+", ");
			}
			out.write(demand.get(demand.size()-1)+"/\n");
			
			out.write("\t"+"dis(pipes) /");
			for(int i=0;i<diameter.size()-1;i++){
				out.write(diameter.get(i)+", ");
			}
			out.write(diameter.get(diameter.size()-1)+"/\n");

			out.write("\t"+"C(pipes) /");
			for(int i=0;i<pipeCost.size()-1;i++){
				out.write(pipeCost.get(i)+", ");
			}
			out.write(pipeCost.get(pipeCost.size()-1)+"/\n");

			out.write("\t"+"R(pipes) /");
			for(int i=0;i<pipeRoughness.size()-1;i++){
				out.write(pipeRoughness.get(i)+", ");
			}
			out.write(pipeRoughness.get(pipeRoughness.size()-1)+"/\n");
			
//			out.write("Source(src)  /");
			
			out.write("\n");
		}
		static void printModel(BufferedWriter out) throws IOException{
			out.write("Scalar omega  /10.68/;\n");
			out.write("Scalar bnd ;\n");
			out.write("Scalar qm;\n");
			out.write("Scalar q_M;\n");

			out.write("\n");

			out.write("bnd = sum(src,D(src));\n");
			out.write("q_M=-bnd;\n");
			out.write("qm=bnd;\n");

			out.write("\n");

			out.write("Variable l(nodes,j,pipes);\n");
			out.write("l.lo(nodes,j,pipes)= 0;\n");

			out.write("\n");

			out.write("Variable q(nodes,j);\n");
			out.write("q.lo(nodes,j)=qm;\n");
			out.write("q.up(nodes,j)=q_M;\n");

			out.write("\n");

			out.write("Variables z;\n");

			out.write("\n");

			out.write("Variable h(nodes);\n");

			out.write("\n");

			out.write("Equations cost \"objective function\",bound1(nodes,j,pipes),cons1(nodes),cons2(nodes),cons3(nodes,j),cons5(src), cons4(nodes,j) ;\n");

			out.write("cost..  z=e=sum(arcs(nodes,j),sum(pipes,l(arcs,pipes)*c(pipes)));\n");

			out.write("bound1(nodes,j,pipes)$arcs(nodes,j).. l(nodes,j,pipes) =l= Len(nodes,j);\n");
			out.write("cons1(nodes).. sum(arcs(j,nodes),q(arcs)) =e= sum(arcs(nodes,j),q(arcs)) + D(nodes);\n");
			out.write("cons2(nodes).. h(nodes) =g= E(nodes) + P(nodes);\n");
			out.write("cons3(arcs(nodes,j)).. h(nodes)-h(j)=e=sum(pipes,((q(arcs)*(abs(q(arcs))**0.852))*(0.001**1.852)*omega*l(arcs,pipes)/((R(pipes)**1.852)*(dis(pipes)/1000)**4.87)));\n");
			out.write("cons4(arcs(nodes,j)).. sum(pipes,l(arcs,pipes)) =e=Len(arcs);\n");
			out.write("cons5(src)..  h(src)=e= sum(srcs,E(srcs));\n");

			out.write("\n");
			
			out.write("model m1  /all/  ;\n");
//			System.out.println("Option threads=4;");
			out.write("m1.optfile =1;\n");
			out.write("solve m1 using minlp minimizing z ;\n");

		}
}
