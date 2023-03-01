package structs;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.Scanner;


public class Testing {
    static final String SOLVER_ROOT_DIR = "/home/manal/Downloads/MTP/JalTantra-Code-and-Scripts";
    // This directory is for temporary use by the method `createNetworkFile()`
    static final String SOLVER_1_NEW_FILE_DIR = "./DataNetworkGraphInput";
    // If `createNetworkFile()` executes successfully, then the created network
    // file will be moved from `SOLVER_1_NEW_FILE_DIR` to this directory
    static final String SOLVER_2_HASH_FILE_DIR = "./DataNetworkGraphInput_hashed";
    // This directory is used by the Python script "CalculateNetworkCost.py"
    // REFER: `OUTPUT_DIR_LEVEL_0` in "CalculateNetworkCost.py"
    static final String SOLVER_3_AUTO_SOLVE_SCRIPT_DIR = "./NetworkResults";

    static final String networkFileHash="525c57589abb519b1dd8e3306fe8dbf5dd11779d1d007854a678db8d7e18171c";
    public static void main(String args[]){
        String Baronm1filePath=SOLVER_ROOT_DIR+"/"+SOLVER_3_AUTO_SOLVE_SCRIPT_DIR+"/"+networkFileHash+"/baron_m1_"+networkFileHash+"/"+"std_out_err.txt";
        String Baronm2filePath=SOLVER_ROOT_DIR+"/"+SOLVER_3_AUTO_SOLVE_SCRIPT_DIR+"/"+networkFileHash+"/baron_m2_"+networkFileHash+"/"+"std_out_err.txt";

        File baronM1File = new File(Baronm1filePath);
        File baronM2File=new File(Baronm2filePath);

        String Remaining_time="";

        String results_for_baron_m1[];
        String cost_for_baron_m1="";

        String results_for_baron_m2[];
        String cost_for_baron_m2="";

        String optimalCost="";

        if(baronM1File.exists()){
            results_for_baron_m1=getIntermediateResults(baronM1File);
            results_for_baron_m2=getIntermediateResults(baronM2File);

            cost_for_baron_m1=results_for_baron_m1[0];
            cost_for_baron_m2=results_for_baron_m2[0];

            Remaining_time=results_for_baron_m1[1];
            if(cost_for_baron_m1.length() > 0){
                System.out.println("baron m1 cost : "+cost_for_baron_m1+" baron m2 cost : "+cost_for_baron_m2 +"time passed : "+Remaining_time);
                if(cost_for_baron_m1.compareTo(cost_for_baron_m2) >= 0){
                    optimalCost=cost_for_baron_m1;
                }
                else optimalCost=cost_for_baron_m2;
            }
        }

    }
    public static String[] getIntermediateResults(File file){
        String res="";
        String Remainingtime="";
        try {
            Scanner scanner = new Scanner(file);
            String lastLine="";
            // read the contents of the file line by line
            while (scanner.hasNextLine()) {
                lastLine = scanner.nextLine();
            }
            lastLine=lastLine.trim();
            System.out.println("last line : "+lastLine);
            String temp[]=lastLine.split("\\s+");
            if(temp.length > 4) {
                try {
                    double number = Double.valueOf(temp[temp.length-1]);
                    double time = Double.valueOf(temp[temp.length-3]);
                    res=""+number;
                    Remainingtime=""+time;
                }
                catch (NumberFormatException e) {
                    System.out.println("There is an exception");
                }
            }
            System.out.println("number : "+res);
            scanner.close();

        } catch (FileNotFoundException e) {
            System.out.println("File not found: " + e.getMessage());
        }
        return new String[]{res,Remainingtime};
    }
}

/***
 *

 boolean resultFileContent=true;
 String resultFilePath=SOLVER_ROOT_DIR + "/" + SOLVER_3_AUTO_SOLVE_SCRIPT_DIR + "/" + networkFileHash + "/0_result.txt";
 String previousRunPath=SOLVER_ROOT_DIR + "/" + SOLVER_3_AUTO_SOLVE_SCRIPT_DIR + "/" + networkFileHash + "/previous_run.txt";

 File resultfile = new File(resultFilePath);
 File previousRunFile=new File(previousRunPath);

 // if somehow the result file was not created
 if (!resultfile.exists()) {
 resultFileContent=false;
 } else {
 // Read from a file
 try (BufferedReader reader = new BufferedReader(new FileReader(resultfile))) {
 String line = reader.readLine();
 if(line.equals("False")) {

 resultFileContent = false;

 // 	delete the previous NetworkResults folder of this particular hash file.
 String pathToFile=SOLVER_ROOT_DIR + "/" + SOLVER_3_AUTO_SOLVE_SCRIPT_DIR + "/" + networkFileHash;
 deleteNetworkResults(pathToFile);

 }
 System.out.println("Result File content: " + line);
 } catch (IOException e) {
 System.out.println("Error reading from Result file: " + e.getMessage());
 }
 System.out.println("Result File already exists.");

 }

 // done to stop the rerunning of that particular network, to avoid the infinte loop, if previously one
 // attempt was made
 if(previousRunFile.exists()){
 // to stop the execution as this has been re run and the results are same as the previous ones
 logd("There was some error in the previous run, delete this network hash and rerun the who;e network");
 System.out.println("This file already existed previously, and there was no error");
 resultFileContent=true;
 }
 else {
 //  create a new file so that next time if the same result comes then
 previousRunFile.createNewFile();
 }
 *
 */
