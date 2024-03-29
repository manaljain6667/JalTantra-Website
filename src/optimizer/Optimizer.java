package optimizer;

import com.google.ortools.linearsolver.MPSolver.ResultStatus;
import optimizer.Pipe.FlowType;
import structs.*;
import javax.servlet.http.HttpServlet;
import javax.servlet.http.HttpServletRequest;
import javax.servlet.http.HttpServletResponse;
import java.io.*;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.security.MessageDigest;
import java.text.ParseException;
import java.text.SimpleDateFormat;
import java.util.*;
import java.util.concurrent.TimeUnit;
import java.util.regex.Matcher;
import java.util.regex.Pattern;




/**
 * This class is responsible for optimization of the network <br> <br>
 */
public class Optimizer {

	/**
	 * Outputs a debug log message.
	 *
	 * @param str the debug log message
	 */
	public static void logd(String str) {
		System.err.println("FM: DEBUG: " + str);
	}

	/**
	 * Outputs an info log message.
	 *
	 * @param str the info log message
	 */
	public static void logi(String str) {
		System.err.println("FM: INFO: " + str);
	}

	/**
	 * Outputs a warning log message.
	 *
	 * @param str the warning log message
	 */
	public static void logw(String str) {
		System.err.println("FM: WARNING: " + str);
	}


	/**
	 * Outputs an error log message.
	 *
	 * @param str the error log message
	 */
	public static void loge(String str) {
		System.err.println("FM: ERROR: " + str);
	}

	public ArrayList<PipeStruct> resultPipes;
	public ArrayList<CommercialPipeStruct> resultCost;
	public ArrayList<ResultPumpStruct> resultPumps;

	private HashMap<Integer, Node> nodes; // nodeID,node
	private HashMap<Integer, Pipe> pipes; // pipeid,pipe
	private List<PipeCost> pipeCost;
	private List<EsrCost> esrCost;

	//following three are data structures received from server request
	private GeneralStruct generalProperties;
	private EsrGeneralStruct esrGeneralProperties;
	private PumpGeneralStruct pumpGeneralProperties;

	private Node source;    //source node of the network

//	private static boolean esrCosting = true;
//	private static boolean esrGen = true;
//	private static boolean removeZeroDemandNodes = false;


	//0 : only pipe optmization, used by default if no esr costing
	//1 : esr costing with only demand nodes allowed
	//2 : general esr with any node allowed
	//3 : gen2 with option to remove certain nodes as esr candidates
	//4 : use l_i_j_k instead of l_i_j
	//5 : better way to compute z_i_j
	//6 : replace how to implement constraint s_i_i = 0 => s_k_i = s_k_p
	//7 : replace above constraint with f_k=1 => sum(s_i_j)=0 and s_i_i=0=>sum(s_i_j)=0
	//8 : added pumps and valves to the optimization
	//9 : pruned some ESR cost rows, depending on the potential downstream demand
	//10: remove all sij variables and instead of node based model, use an edge based model
	private int modelNumber = 0;


	//default values are provided for following parameters, but are typically overridden by user in the server request

	// flow in secondary = factor * flow in primary
	// primary pumping hours = factor * secondary pumping hours
	private double secondaryFlowFactor = 2;

	// ratio of size of ESR to the daily demand of the nodes it serves
	private double esrCapacityFactor = 1;

	// minimum and maximum height allowed for an ESR
	private double maxEsrHeight = 25;

	// minimum and maximum power allowed for a pump
	private double minPumpPower = 1;
	private double maxPumpPower = 10000;

	// maximum pressure head that can be provided by a pump
	private int maxPumpHead = 10000;

	//container for pump and valve information
	private PumpManualStruct[] pumpManualArray;
	private ValveStruct[] valves;

	//create an instance of optimizer
	/**
	 * Constructs an instance of the {@link Optimizer} class.
	 *
	 * @param nodeStructs            an array of {@link NodeStruct} objects representing the node properties
	 * @param pipeStructs            an array of {@link PipeStruct} objects representing the pipe properties
	 * @param commercialPipeStructs  an array of {@link CommercialPipeStruct} objects representing the commercial pipe properties
	 * @param generalStruct          a {@link GeneralStruct} object representing the general network properties
	 * @param esrGeneralProperties   an {@link EsrGeneralStruct} object representing the ESR general properties (can be null if ESR optimization is disabled)
	 * @param esrCostsArray          an array of {@link EsrCostStruct} objects representing the ESR cost properties (can be null if ESR optimization is disabled)
	 * @param pumpGeneralProperties  a {@link PumpGeneralStruct} object representing the pump general properties (can be null if pump optimization is disabled)
	 * @param pumpManualArray        an array of {@link PumpManualStruct} objects representing the manual pump properties (can be null if pump optimization is disabled)
	 * @param valves                 an array of {@link ValveStruct} objects representing the valve properties
	 * @throws Exception if there is an error in the construction of the optimizer
	 */
	public Optimizer(NodeStruct[] nodeStructs, PipeStruct[] pipeStructs, CommercialPipeStruct[] commercialPipeStructs, GeneralStruct generalStruct, EsrGeneralStruct esrGeneralProperties, EsrCostStruct[] esrCostsArray, PumpGeneralStruct pumpGeneralProperties, PumpManualStruct[] pumpManualArray, ValveStruct[] valves) throws Exception {
		nodes = new HashMap<Integer, Node>();
		pipes = new HashMap<Integer, Pipe>();
		pipeCost = new ArrayList<PipeCost>();
		esrCost = new ArrayList<EsrCost>();
		generalProperties = generalStruct;

		int[] a = {};
		this.pumpGeneralProperties = new PumpGeneralStruct(false, 1, 100, 0, 0, 1, 0, 0, 1, a);

		Set<Integer> usedNodeIDs = new HashSet<Integer>();

		// Initialize source node
		// NOTE: Pressure at each node = "Head - Elevation"
		//       However, minimum pressure required at the source is always set to 0 irrespective of the
		//       value of `generalProperties.min_node_pressure`. This is based on network file analysis.
		//       For more details, have a look at `createNetworkFile()` method where "E = Elevation of
		//       each node" is being written to the output file.
		source = new Node(
				generalProperties.source_elevation,
				0.0,
				generalProperties.source_nodeid,
				0,  // generalProperties.source_head - generalProperties.source_elevation
				generalProperties.source_nodename,
				24 / generalProperties.supply_hours,
				usedNodeIDs
		);
		source.setAllowESR(false);
		usedNodeIDs.add(source.getNodeID());
		source.setHead(generalProperties.source_head);
		nodes.put(source.getNodeID(), source);

		//initialize all the other nodes
		double totalDemand = 0.0;
		for (NodeStruct node : nodeStructs) {
			double minPressure = node.minpressure == 0 ? generalProperties.min_node_pressure : node.minpressure;
			Node n = new Node(node.elevation, node.demand, node.nodeid, minPressure, node.nodename, 24 / generalProperties.supply_hours, usedNodeIDs);
			usedNodeIDs.add(n.getNodeID());
			nodes.put(n.getNodeID(), n);
			totalDemand += n.getDemand();
		}
		// NOTE: Demand of source node = -1 * sum(demand of other nodes)
		//       This is based on network file analysis
		// nodes.get(source.getNodeID()).setDemand(-totalDemand);
		this.source.setDemand(-totalDemand);  // This will update the element pointing to the source node in `nodes` HashMap as well

		Set<Integer> usedPipeIDs = new HashSet<Integer>();

		//initialize the pipes
		for (PipeStruct pipe : pipeStructs) {
			double roughness = pipe.roughness == 0 ? generalProperties.def_pipe_roughness : pipe.roughness;

			Node startNode = nodes.get(pipe.startnode);
			if (startNode == null) {
				throw new Exception("Invalid startNode:" + pipe.startnode + " provided for pipe ID:" + pipe.pipeid);
			}

			Node endNode = nodes.get(pipe.endnode);
			if (endNode == null) {
				throw new Exception("Invalid endNode:" + pipe.endnode + " provided for pipe ID:" + pipe.pipeid);
			}

			Pipe p = new Pipe(pipe.length, startNode, endNode, pipe.diameter, roughness, pipe.pipeid, pipe.parallelallowed, usedPipeIDs);
			usedPipeIDs.add(p.getPipeID());
			pipes.put(p.getPipeID(), p);
		}

		//initialize the commercial pipe information
		for (CommercialPipeStruct commercialPipe : commercialPipeStructs) {
			double roughness = commercialPipe.roughness == 0 ? generalProperties.def_pipe_roughness : commercialPipe.roughness;
			pipeCost.add(new PipeCost(commercialPipe.diameter, commercialPipe.cost, Double.MAX_VALUE, roughness));
		}
		this.valves = valves;

		//default model number is 0 for only pipe optimization
		modelNumber = 0;

		//if ESR optimization enabled, initialize ESR properties and set modelnumber
		if (esrGeneralProperties != null && esrGeneralProperties.esr_enabled) {
			this.esrGeneralProperties = esrGeneralProperties;

			if (esrGeneralProperties.secondary_supply_hours == 0) {
				throw new Exception("ESR option is enabled, but secondary supply hours is provided as zero.");
			}

			if (esrGeneralProperties.esr_capacity_factor == 0) {
				throw new Exception("ESR option is enabled, but esr capacity factor is provided as zero.");
			}

			secondaryFlowFactor = generalProperties.supply_hours / esrGeneralProperties.secondary_supply_hours;
			esrCapacityFactor = esrGeneralProperties.esr_capacity_factor;
			maxEsrHeight = esrGeneralProperties.max_esr_height;

			modelNumber = 9;

			for (EsrCostStruct esrcost : esrCostsArray) {
				esrCost.add(new EsrCost(esrcost.mincapacity,
						esrcost.maxcapacity,
						esrcost.basecost,
						esrcost.unitcost));
			}
		}

		//if pump enabled, initialize pump properties
		if (pumpGeneralProperties != null && pumpGeneralProperties.pump_enabled) {
			this.pumpGeneralProperties = pumpGeneralProperties;
			this.pumpManualArray = pumpManualArray;
			this.minPumpPower = pumpGeneralProperties.minpumpsize;

			if (pumpGeneralProperties.efficiency == 0)
				throw new Exception("Pump option is enabled, but pump efficiency is provided as zero.");

			if (pumpGeneralProperties.design_lifetime == 0)
				throw new Exception("Pump option is enabled, but design lifetime is provided as zero.");
		}

		//set total demand required for the network
		totalDemand = getTotalCapacity();
	}

	/**
	 * Returns the general properties of the optimizer.
	 *
	 * @return the general properties of the optimizer
	 */
	public GeneralStruct getGeneralProperties() {
		return generalProperties;
	}

	/**
	 * Returns the ESR (Elevated Storage Reservoir) general properties of the optimizer.
	 *
	 * @return the ESR general properties of the optimizer
	 */
	public EsrGeneralStruct getEsrGeneralProperties() {
		return esrGeneralProperties;
	}

	/**
	 * Returns the pump general properties of the optimizer.
	 *
	 * @return the pump general properties of the optimizer
	 */
	public PumpGeneralStruct getPumpGeneralProperties() {
		return pumpGeneralProperties;
	}

	/**
	 * Returns the project name from the general properties of the optimizer.
	 *
	 * @return the project name from the general properties
	 */
	public String getGeneralPropertiesProjectName() {
		return generalProperties.name_project;
	}

	/**
	 * Returns the organization name from the general properties of the optimizer.
	 *
	 * @return the organization name from the general properties
	 */
	public String getGeneralPropertiesOrganizationName() {
		return generalProperties.name_organization;
	}


	/**
	 * Check whether the network structure is valid or not
	 *
	 * @return int <br>
	 * 1 => Valid Input: Acyclic graph<br>
	 * 2 => Valid Input: Cyclic graph<br>
	 * 3 => Invalid Input: Disconnected graph (i.e. nodes disconnected from the network)<br>
	 * 4 => Invalid Input: Source Head is less than Source Elevation
	 */
	private int validateNetwork() {
		if (!(this.source.getHead() >= this.source.getElevation())) {
			return 4;
		}

		Node root = source;
		HashSet<Node> seen = new HashSet<>();
		Stack<Node> left = new Stack<>();
		left.add(root);

		while (!left.isEmpty()) {
			Node top = left.pop();
			if (seen.contains(top)) {
				return 2; // cycle
			}
			seen.add(top);
			for (Pipe pipe : top.getOutgoingPipes()) {
				left.push(pipe.getEndNode());
			}
		}

		if (seen.size() != nodes.size()) {
			System.out.println(seen.size());
			System.out.println(nodes.size());
			return 3;  // not fully connected
		}
		return 2;
	}

	//return the total ESR capacity required in the network in litres

	/**
	 * Calculates and returns the total capacity required for the network.
	 *
	 * @return the total capacity required for the network
	 */
	private double getTotalCapacity() {
		double sum = 0;
		for (Node n : nodes.values()) {
			sum = sum + n.getRequiredCapacity(esrCapacityFactor);
		}
		return sum;
	}

	// ---

	public static int attempts=0;

	// static final String SOLVER_ROOT_DIR = "/home/deploy";

	/**
	 * This should point to the root of the repository: Jaltantra-Code-and-Scripts
	 */
	static final String SOLVER_ROOT_DIR = "/home/manal/Downloads/MTP/JalTantra-Code-and-Scripts";
	// This directory is for temporary use by the method `createNetworkFile()`
	static final String SOLVER_1_NEW_FILE_DIR = "./DataNetworkGraphInput";
	// If `createNetworkFile()` executes successfully, then the created network
	// file will be moved from `SOLVER_1_NEW_FILE_DIR` to this directory
	static final String SOLVER_2_HASH_FILE_DIR = "./DataNetworkGraphInput_hashed";
	// This directory is used by the Python script "CalculateNetworkCost.py"
	// REFER: `OUTPUT_DIR_LEVEL_0` in "CalculateNetworkCost.py"
	static final String SOLVER_3_AUTO_SOLVE_SCRIPT_DIR = "./NetworkResults";

	// Amount of time "CalculateNetworkCost.py" should execute the solver(s) for the network file
	static String SOLVER_EXECUTION_TIME = "00:05:00";
	static String SOLVER_EXECUTION_TIME_DISPLAY_STR = "5 minutes";

	static String list_of_times[]={"00:10:00","00:20:00"};

	static String run_Time="";


	/**
	 * Optimize the network. And, if done successfully, the results are stored in three ArrayLists,
	 * namely resultPipes, resultCost and resultPumps.
	 * @param runTime
	 * 			a string variable denoting the run time of the network
	 * @return whether network was solved successfully results are ready (`true`) or not (`false`)
	 * @throws Exception in case of any error/problem or to convey any status information
	 */
	public boolean Optimize(String runTime) throws Exception {

		// Execution flow:
		//   1. Create the data files for the network
		//   2. Asynchronously launch `CalculateNetworkCost.py` for the data files

		run_Time=runTime;

		if(runTime.equals("1hour")){
			SOLVER_EXECUTION_TIME="01:00:00";
			SOLVER_EXECUTION_TIME_DISPLAY_STR="1 hour";
		}

		if(runTime.equals("5min")){
			SOLVER_EXECUTION_TIME="00:05:00";
			SOLVER_EXECUTION_TIME_DISPLAY_STR="5 minutes";
		}

		// Validate the network layout
		logd("Network validation started...");
		final int networkValidationResult = validateNetwork();
		logd("Network validation complete..., networkValidationResult = " + networkValidationResult);
		if (networkValidationResult == 1 || networkValidationResult == 2) {
			final String networkFileResult = createNetworkFile();
			final String networkFileStatus = networkFileResult.substring(0, networkFileResult.indexOf("-"));
			final String networkFileName = networkFileResult.substring(2);
			String networkFileHash = networkFileName.substring(0, networkFileName.lastIndexOf("."));

			/**
			 *
			 here prompt the user to save the network file before clicking the optimize button
			 logic:
			 - create a temp file to see whether this is the first time to prompt the user,
			 if it exists then do not prompt the user otherwise prompt the user
			 *
			 */
			String promptFilePAth=SOLVER_ROOT_DIR+"/"+SOLVER_2_HASH_FILE_DIR+"/"+networkFileHash+"_prompt.txt";
			boolean promptFileExist=checkFileExist(promptFilePAth);

			if(!promptFileExist){
				logi("Please save the network file, if not, before continuing the optimization");
				throw new Exception("Please save the network file, if not, before continuing the optimization and click on optimize button again");
			}

			String previousHash=networkFileHash;

			networkFileHash=networkFileHash+runTime;

			boolean statusFileExists = (new File(SOLVER_ROOT_DIR + "/" + SOLVER_3_AUTO_SOLVE_SCRIPT_DIR + "/" + networkFileHash + "/0_status")).exists();

			// if the statusFileExists, but previously there was no output and the content of this file was false
			// then we need to re run the jaltantra for this network
			boolean resultFileContent=true;

			String resultFilePath=SOLVER_ROOT_DIR + "/" + SOLVER_3_AUTO_SOLVE_SCRIPT_DIR + "/" + networkFileHash + "/0_result.txt";
			String previousRunPath=SOLVER_ROOT_DIR + "/" + SOLVER_2_HASH_FILE_DIR + "/" + networkFileName + "_previous_run.txt";
			int cntlinesInResultFile=0;

			File resultfile = new File(resultFilePath);
			File previousRunFile=new File(previousRunPath);

			String message="";

			logd("networkFileHash = " + networkFileHash + ", statusFileExists = " + statusFileExists);

			if (networkFileStatus.equals("0") || statusFileExists == false ) {
				logi("Starting the solvers (CalculateNetworkCost.py) for network file '" + networkFileResult + "'");
				launchCalculateNetworkCost(SOLVER_ROOT_DIR + "/" + SOLVER_2_HASH_FILE_DIR + "/" + networkFileName);
				throw new Exception("The solver is working, refresh the page after " + SOLVER_EXECUTION_TIME_DISPLAY_STR + " to see the results");
			} else if (networkFileStatus.equals("2")) {
				// Check which of the following case is true:
				//   - solvers are running
				//   - solvers are done executing, but with an error
				//   - results already generated
				int status = checkSolverResultStatus(SOLVER_ROOT_DIR + "/" + SOLVER_3_AUTO_SOLVE_SCRIPT_DIR + "/" + networkFileHash);
				if (status == -1) {
					loge("CHECKME: FIXME: Probably failed to start the solver launch script, or 'CalculateNetworkCost.py' failed to create the status file");

					// here we will check for the launching error of the file, if there was then delete the previous directory

					String pathToFile=SOLVER_ROOT_DIR + "/" + SOLVER_3_AUTO_SOLVE_SCRIPT_DIR + "/" + networkFileHash;
					deleteNetworkResults(pathToFile);

					throw new Exception("Internal server error: failed to start the solver launch script, click on optimize button again");
				} else if (status == 0) {
					loge("FIXME: Failed to start the solver for network file: " + networkFileResult);
					throw new Exception("Internal server error: failed to start the solver for this network");
				} else if (status == 1) {
					String Baronm1filePath=SOLVER_ROOT_DIR+"/"+SOLVER_3_AUTO_SOLVE_SCRIPT_DIR+"/"+networkFileHash+"/baron_m1_"+previousHash+"/"+"std_out_err.txt";
					String Baronm2filePath=SOLVER_ROOT_DIR+"/"+SOLVER_3_AUTO_SOLVE_SCRIPT_DIR+"/"+networkFileHash+"/baron_m2_"+previousHash+"/"+"std_out_err.txt";

					File baronM1File = new File(Baronm1filePath);
					File baronM2File=new File(Baronm2filePath);

					String time_passed="";

					String results_for_baron_m1[];
					String cost_for_baron_m1="";

					String results_for_baron_m2[];
					String cost_for_baron_m2="";

					String optimalCost="";

					results_for_baron_m1=getIntermediateResults(baronM1File);
					results_for_baron_m2=getIntermediateResults(baronM2File);

					cost_for_baron_m1=results_for_baron_m1[0];
					cost_for_baron_m2=results_for_baron_m2[0];

					// done as baron m2 is processed later
					time_passed=results_for_baron_m2[1];
					if(time_passed.length() == 0) time_passed=results_for_baron_m1[1];
					if(cost_for_baron_m1.length() > 0){
						System.out.println("baron m1 cost : "+cost_for_baron_m1+" baron m2 cost : "+cost_for_baron_m2 +"time passed : "+time_passed);
						if(cost_for_baron_m1.compareTo(cost_for_baron_m2) >= 0){
							optimalCost=cost_for_baron_m2;
						}
						else optimalCost=cost_for_baron_m1;
					}

					String log_file=SOLVER_ROOT_DIR +"/" + SOLVER_2_HASH_FILE_DIR + "/" + previousHash + ".R" + runTime + ".log";
					File log_file_path = new File(log_file);
					String time_to_display[] = getIntermediateTime(log_file_path);

					if(time_to_display[0].length() > 0) {

						String timeFinished = time_to_display[0];
						String timeLeft = time_to_display[1];
						logi("Solver is running for this network. Please wait..., current cost : "+optimalCost);
						throw new Exception("Solver is running for this network. Please wait...,  "+timeFinished+" passed, remaining Time : "+timeLeft);

						// convert sec into hours minutes and seconds format
						// String resultTime=getTimeToDisplay(time_passed);
						//
						// if(runTime.equals("5min")){
						//
						// 	String remainingTime="";
						// 	String pathToattemptsFile=SOLVER_ROOT_DIR+"/"+SOLVER_2_HASH_FILE_DIR+"/"+networkFileName+"_attemptsFile.txt";
						//
						// 	if(!new File(pathToattemptsFile).exists()){
						// 		int remaining=300-((int)(Float.parseFloat(time_passed)));
						// 		remainingTime=getTimeToDisplay(""+remaining);
						// 	}
						// 	else{
						// 		int temp=0;
						// 		try (BufferedReader reader = new BufferedReader(new FileReader(pathToattemptsFile))) {
						// 			String line = reader.readLine();
						// 			temp=Integer.parseInt(line);
						// 		} catch (IOException e) {
						// 			System.out.println("Error reading from attempts file during displaying of time: " + e.getMessage());
						// 		}
						// 		String currentTime=list_of_times[temp-1];
						// 		String currentMinutes=getMinutesToDisplay(currentTime);
						//
						// 		int currentMin=Integer.parseInt(currentMinutes);
						// 		int currentsec=currentMin*60;  // we will get seconds
						//
						// 		int remining=currentsec-((int)(Float.parseFloat(time_passed)));
						// 		remainingTime=getTimeToDisplay(""+remining);
						// 	}
						// 	logi("Solver is running for this network. Please wait..., current cost : "+optimalCost);
						// 	throw new Exception("Solver is running for this network. Please wait...,  "+resultTime+" passed, remaining Time : "+remainingTime);
						// }
						// logi("Solver is running for this network. Please wait..., current cost : "+optimalCost);
						// throw new Exception("Solver is running for this network. Please wait...,  "+resultTime+" passed");
					}
					else{
						logi("Solver is running for this network. Please wait...");
						throw new Exception("Solver is running for this network. Please wait...");
					}

					// logi("Solver is running for this network. Please wait...");
					// throw new Exception("Solver is running for this network. Please wait...");
				} else if (status == 2) {

					if(!runTime.equals("5min")){
						throw new Exception("Either no feasible solution found, or failed to solve the network in " + SOLVER_EXECUTION_TIME_DISPLAY_STR+" ");
					}

					loge("CHECKME: The solvers finished the execution, but failed to get the result. " +
						 "Either some unknown error, or no feasible solution found, or failed to solve the network, click optimize again and check after " + SOLVER_EXECUTION_TIME_DISPLAY_STR+" minutes");
					//creating an attempts file and relaunching the network with the updated time
					String pathToattemptsFile=SOLVER_ROOT_DIR+"/"+SOLVER_2_HASH_FILE_DIR+"/"+networkFileName;
					File file = new File(pathToattemptsFile+"_attemptsFile.txt");
					if (!file.exists()) {
						try {
							file.createNewFile();
							// Write to a file
							String content = "0";
							attempts=0;
							try (FileWriter writer = new FileWriter(file)) {
								writer.write(content);
								System.out.println("File written successfully.");
							} catch (IOException e) {
								System.out.println("Error writing to file: " + e.getMessage());
							}
							System.out.println("File created successfully.");
						} catch (IOException e) {
							System.out.println("Error creating file: " + e.getMessage());
						}
					} else {
						// Read from a file
						try (BufferedReader reader = new BufferedReader(new FileReader(file))) {
							String line = reader.readLine();
							attempts=Integer.parseInt(line);
							System.out.println("File content: " + line);
						} catch (IOException e) {
							System.out.println("Error reading from file: " + e.getMessage());
						}
						System.out.println("File already exists.");
					}

					for(int i=attempts;i<list_of_times.length;i++){
						SOLVER_EXECUTION_TIME=list_of_times[attempts];
						attempts++;
						try (FileWriter writer = new FileWriter(file)) {
							writer.write(""+attempts);
							System.out.println("File written successfully.");
						} catch (IOException e) {
							System.out.println("Error writing to file: " + e.getMessage());
						}
						System.out.println("system is running for "+SOLVER_EXECUTION_TIME+"minutes, attempt number "+attempts);

						//delete previously existed network results
						String pathToFile=SOLVER_ROOT_DIR + "/" + SOLVER_3_AUTO_SOLVE_SCRIPT_DIR + "/" + networkFileHash;
						deleteNetworkResults(pathToFile);
						statusFileExists = (new File(SOLVER_ROOT_DIR + "/" + SOLVER_3_AUTO_SOLVE_SCRIPT_DIR + "/" + networkFileHash + "/0_status")).exists();
						launchCalculateNetworkCost(SOLVER_ROOT_DIR + "/" + SOLVER_2_HASH_FILE_DIR + "/" + networkFileName);
						status = checkSolverResultStatus(SOLVER_ROOT_DIR + "/" + SOLVER_3_AUTO_SOLVE_SCRIPT_DIR + "/" + networkFileHash);
						// System.out.println("system is now executing for 10 minutes");
						if (status == 3) {
							logi("Extracting the result");
							boolean ok = extractSolverResult(SOLVER_ROOT_DIR + "/" + SOLVER_3_AUTO_SOLVE_SCRIPT_DIR + "/" + networkFileHash);
							if (ok) return true;
							loge("CHECKME: FIXME: extractSolverResult(...) return false for network file with hash =" + networkFileHash);
							String pathToFolder=SOLVER_ROOT_DIR + "/" + SOLVER_3_AUTO_SOLVE_SCRIPT_DIR + "/" + networkFileHash;
							deleteNetworkResults(pathToFolder);
							throw new Exception("Internal server error: result extraction failed for the network file with hash =" + networkFileHash+" click optimize button again or after some time to rerun the file");
						}
						SOLVER_EXECUTION_TIME_DISPLAY_STR=""+getMinutesToDisplay(SOLVER_EXECUTION_TIME);
						throw new Exception("Either no feasible solution found, or failed to solve the network, retrying for " + SOLVER_EXECUTION_TIME_DISPLAY_STR+" minutes");
					}
					if(attempts > 0){
						SOLVER_EXECUTION_TIME_DISPLAY_STR=""+getMinutesToDisplay(list_of_times[attempts-1]);
					}
					throw new Exception("Either no feasible solution found, or failed to solve the network in " + SOLVER_EXECUTION_TIME_DISPLAY_STR+" minutes");

				} else if (status == 3) {
					logi("Extracting the result");
					boolean ok = extractSolverResult(SOLVER_ROOT_DIR + "/" + SOLVER_3_AUTO_SOLVE_SCRIPT_DIR + "/" + networkFileHash);
					if (ok) return true;
					loge("CHECKME: FIXME: extractSolverResult(...) return false for network file with hash =" + networkFileHash);
					String pathToFolder=SOLVER_ROOT_DIR + "/" + SOLVER_3_AUTO_SOLVE_SCRIPT_DIR + "/" + networkFileHash;
					deleteNetworkResults(pathToFolder);
					throw new Exception("Internal server error: result extraction failed for the network file with hash =" + networkFileHash+" click optimize button again or after some time to get the results");
				} else {
					loge("CHECKME: FIXME: Unexpected return value from `checkSolverResultStatus()`. deleted the past execution result for hash =" + networkFileHash+" , click the optimize button again");

					String pathToFolder=SOLVER_ROOT_DIR + "/" + SOLVER_3_AUTO_SOLVE_SCRIPT_DIR + "/" + networkFileHash;
					deleteNetworkResults(pathToFolder);

					throw new Exception("Internal server error: unknown execution status. Need to delete the past execution result for hash =" + networkFileHash+" , click the optimize button again");
				}
			} else {
				loge("FIXME: `createNetworkFile()` method returned unexpected value ' " + networkFileResult + " '");
				throw new Exception("Internal server error: createNetworkFile() method failed, probably due to unexpected change in some method(s)");
			}
		}
		else if (networkValidationResult == 3) {
			throw new Exception("Input is not valid. Nodes unconnected in the network");
		} else if (networkValidationResult == 4) {
			throw new Exception("Source Head (" + this.source.getHead() + ") should be greater than or equal to Source Elevation (" + this.source.getElevation() + "). Please fix it in the 'General' section");
		} else {
			loge("FIXME: network validation returned unexpected value " + networkValidationResult);
			throw new Exception("Internal server error: network validation failed");
		}
		// There is NO case in which the program will reach this point
		// return true;
	}

	/**
	 * This function takes a string value which denotes seconds and converts it to HH:MM:SS string format
	 * @param timeInSec
	 * 			a string variable denoting the time in seconds
	 * @return converted value of time from seconds to HH:MM:SS format
	 */
	public static String getTimeToDisplay(String timeInSec){
		String resultTime="";
		int sec=(int)(Float.parseFloat(timeInSec));

		int hours=sec/3600;
		sec=sec%3600;

		int minutes=sec/60;
		sec=sec%60;

		if(hours > 0) resultTime=""+hours+"hours ";
		if(minutes > 0) resultTime=""+minutes+"minutes "+sec+" seconds";
		else resultTime="00 minutes"+sec+" seconds";

		return resultTime;
	}

	/**
	 * This function checks whether a particular file exist or not. If not exist then it creates that file<br>
	 * in that particular filePath
	 * @param filePath
	 * 		a string variable which denotes the path to the file
	 * @return boolean <br>
	 * 	true => if the file exists
	 * 	false => if that file does not exist and created that file
	 */
	public static boolean checkFileExist(String filePath){

		boolean exists=false;

		File file = new File(filePath);
		if (!file.exists()) {
			try {
				file.createNewFile();
			} catch (IOException e) {
				System.out.println("Error creating file: " + e.getMessage());
			}
		}
		else exists=true;

		return exists;
	}

	/**
	 * Retrieves the intermediate time information from the given file.
	 *
	 * @param file The file to read and extract the intermediate time information from.
	 * @return An array of two strings containing the time finished and time left respectively.
	 */
	public static String[] getIntermediateTime(File file){
		String timeFinished="";
		String timeLeft="";
		try {
			Scanner scanner = new Scanner(file);
			String lastLine="";
			// read the contents of the file line by line
			while (scanner.hasNextLine()) {
				String line = scanner.nextLine();
				Pattern pattern = Pattern.compile("Time Finished =");
				Matcher matcher = pattern.matcher(line);
				if(matcher.find()){
					lastLine = line;
				}
			}
			lastLine=lastLine.trim();
			System.out.println("last line : "+lastLine);

			Pattern pattern = Pattern.compile("Time Finished = (\\d{1,2}:\\d{2}:\\d{2}), Time Left = (\\d{1,2}:\\d{2}:\\d{2})");
			Matcher matcher = pattern.matcher(lastLine);
			if (matcher.find()) {
				timeFinished = matcher.group(1);
				timeLeft = matcher.group(2);
				System.out.println("Time Finished: " + timeFinished);
				System.out.println("Time Left: " + timeLeft);
			}
			System.out.println("Time finished : " + timeFinished + "Time left : " + timeLeft);
			scanner.close();

		} catch (FileNotFoundException e) {
			System.out.println("Either File not found or not able to extract time: " + e.getMessage());
		}
		return new String[]{timeFinished , timeLeft};
	}

	/**
	 * Retrieves the intermediate results from the given file.
	 *
	 * @param file The file to read and extract the intermediate results from.
	 * @return An array of two strings containing the result value and remaining time respectively.
	 */
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

	/**
	 * Deletes the results of a network stored at the specified path.
	 *
	 * @param pathToHashedNetwork The path to the hashed network directory.
	 */
	public void deleteNetworkResults(String pathToHashedNetwork){
		File directory = new File(pathToHashedNetwork);

		if (directory.isDirectory()) {
			File[] files = directory.listFiles();
			for (File file : files) {
				file.delete();
			}
			directory.delete();
			System.out.println("Directory deleted successfully.");
		} else {
			System.out.println("Error: Not a directory.");
		}
	}

	// public static String parseDate(String dateString,int multiplier) {
	// 	SimpleDateFormat dateFormat = new SimpleDateFormat("HH:mm:ss");
	// 	String updatedTime="";
	// 	try {
	// 		System.out.println(dateString);
	// 		Date date = dateFormat.parse(dateString);
	//
	// 		int doubled=date.getMinutes() * multiplier;
	// 		System.out.println("minutes : "+ doubled);
	//
	// 		int hours = doubled / 60;
	// 		int remainingMinutes = doubled % 60;
	// 		int seconds = 0;
	//
	// 		updatedTime = String.format("%02d:%02d:%02d", hours, remainingMinutes, seconds);
	// 		System.out.println("Time: " + updatedTime);
	//
	// 	} catch (ParseException e) {
	// 		System.out.println("Error parsing date string: " + e.getMessage());
	// 	}
	//
	// 	return updatedTime;
	// }

	/**
	 *
	 * Converts a given date string into minutes and returns it as a string.
	 *
	 * @param dateString The date string in the format "HH:mm:ss".
	 * @return The minutes extracted from the date string as a string.
	 */
	public String getMinutesToDisplay(String dateString) {
		SimpleDateFormat dateFormat = new SimpleDateFormat("HH:mm:ss");
		int minutes=0;
		try {
			System.out.println(dateString);
			Date date = dateFormat.parse(dateString);

			minutes=date.getMinutes();
			System.out.println("minutes : "+ minutes);


		} catch (ParseException e) {
			System.out.println("Error parsing date string: " + e.getMessage());
		}

		return ""+minutes;
	}

	/**
	 * Create directory if it does not exist.
	 *
	 * @param dirPath path to the directory which is to be created
	 * @throws Exception if directory does not exist and its creation failed
	 */
	private static void checkAndCreateDir(final String dirPath) throws Exception {
		// REFER: https://stackoverflow.com/questions/3634853/how-to-create-a-directory-in-java
		File theDir = new File(dirPath);
		if (!theDir.exists()) {
			boolean res = theDir.mkdirs();
			if (!res) {
				loge("FIXME: Failed to create directory: '" + dirPath + "'");
				throw new Exception("Internal error at the server, the backend does not have write permission");
			}
		}
	}

	/**
	 * Generate unique filename based on the time at which the function was called and a
	 * random number, such that no file with the same name exists in `baseDirectoryPath`.
	 * <br>
	 * The filename format is: "yyyy-MM-dd_HH-mm-ss.SSS_{i}_{RandomNumber}.R"
	 * Range of i = [0, 9_99_999]
	 * Range of RandomNumber = [10_00_001, 99_99_999]
	 *
	 * @param baseDirectoryPath path to the directory for which a unique filename is to be generated
	 * @return unique filename which does not exist
	 * @throws Exception if unique filename generation fails
	 */
	private static String generateUniqueFileName(final String baseDirectoryPath) throws Exception {
		checkAndCreateDir(baseDirectoryPath);

		// REFER: https://www.java-examples.com/formatting-date-custom-formats-using-simpledateformat
		// REFER: https://docs.oracle.com/javase/7/docs/api/java/text/SimpleDateFormat.html
		Date date = new Date();
		SimpleDateFormat sdf = new SimpleDateFormat("yyyy-MM-dd_HH-mm-ss.SSS_");
		String strDate = sdf.format(date);

		// REFER: https://stackoverflow.com/questions/5887709/getting-random-numbers-in-java
		// 7 digit random number with range [1000001, 9999999]
		long randomValue = (long) (1000000 + (8999999 * Math.random() + 1));

		for (int i = 0; i < 1000000; ++i) {
			if ((new File(baseDirectoryPath + "/" + strDate + i + "_" + randomValue + ".R")).exists())
				continue;
			return strDate + i + "_" + randomValue + ".R";
		}

		// REFER: https://stackoverflow.com/questions/840190/changing-the-current-working-directory-in-java
		// REFER: https://stackoverflow.com/questions/4871051/how-to-get-the-current-working-directory-in-java
		logi("pwd: " + System.getProperty("user.dir"));
		loge("CHECKME: FIXME: Could not generate unique filename, date time = " + strDate);
		throw new Exception(
				"Internal error at the server, could not generate a unique filename for the " +
				"request (at " + strDate + "). Probably the server is too loaded. Please try again."
		);
	}

	/**
	 * This method return the digest/hash (in hexadecimal format) of the file passed<br>
	 * REFER: https://www.geeksforgeeks.org/how-to-generate-md5-checksum-for-files-in-java/
	 *
	 * @param digest hashing algorithm to be used
	 * @param file   file object whose hash is to be found
	 * @return Hex representation of the hash
	 * @throws IOException if file read operation fails
	 */
	private static String findFileHashInHex(MessageDigest digest, File file) throws IOException {
		// Get file input stream for reading the file content
		FileInputStream fis = new FileInputStream(file);

		// Create byte array to read data in chunks
		byte[] byteArray = new byte[1024];
		int bytesCount;

		// Read the data from file and update that data in the message digest
		while ((bytesCount = fis.read(byteArray)) != -1) {
			digest.update(byteArray, 0, bytesCount);
		}

		// Close the input stream
		fis.close();

		// Store the bytes returned by the digest() method
		byte[] bytes = digest.digest();

		// This array of bytes has bytes in decimal format, so we need to convert it into hexadecimal format
		// For this we create an object of StringBuilder since it allows us to update the string i.e. its mutable
		StringBuilder sb = new StringBuilder();

		// Loop through the bytes array
		for (byte aByte : bytes) {
			// The following line converts the decimal into hexadecimal format and appends that to the
			// StringBuilder object
			sb.append(Integer.toString((aByte & 0xff) + 0x100, 16).substring(1));
		}



		// Finally, we return the complete hash
		return sb.toString().toLowerCase();
	}

	/**
	 * Read a file completely and return its content.
	 *
	 * @param pathToFile path to the file which is to be read
	 * @return content of the file as String
	 * @throws IOException if readAllBytes(...) fails
	 */
	private static String readFileAsString(final String pathToFile) throws IOException {
		// REFER: https://www.geeksforgeeks.org/different-ways-reading-text-file-java/
		return new String(Files.readAllBytes(Paths.get(pathToFile)));
	}

	/**
	 * Returns the MD5 hash of the newly created network file provided everything goes properly
	 *
	 * @return String<br>
	 * "0-NetworkFileMD5Hash.R" if everything was fine and this is the first request for this network<br>
	 * "2-NetworkFileMD5Hash.R" if everything was fine and a request for this network was received in the past
	 * @throws Exception if network file creation fails due to any reason
	 */
	private String createNetworkFile() throws Exception {
		final String uniqueFileName = generateUniqueFileName(SOLVER_ROOT_DIR + "/" + SOLVER_1_NEW_FILE_DIR);
		final String pathTo1FreshFile = SOLVER_ROOT_DIR + "/" + SOLVER_1_NEW_FILE_DIR + "/" + uniqueFileName;

		// REFER: https://stackoverflow.com/questions/10667734/java-file-open-a-file-and-write-to-it
		BufferedWriter out = null;
		try {
			// Network file creation failed as a file already exists with the same name.
			// Probably, the server is loaded with too many requests.
			if ((new File(pathTo1FreshFile)).exists())
				throw new Exception("Network file creation failed because the server has too many requests running. Please retry after a few moments.");

			out = new BufferedWriter(new FileWriter(pathTo1FreshFile));

			// TODO: NOT SURE, but, this may need to be updated if we have to use `nodeName` instead of `nodeID`
			// NOTE: Set of nodes/vertexes
			out.write("set nodes :=");
			for (Node n : this.nodes.values()) {
				out.write(" " + n.getNodeID());
			}
			out.write(";");

			// NOTE: Set of commercial pipes available
			out.write("\n\nset pipes :=");
			for (int i = 0; i < this.pipeCost.size(); ++i) {
				out.write(" " + i);
			}
			out.write(";");

			// NOTE: Set of arcs/links/edges
			out.write("\n\nparam : arcs : L :=");  // NOTE: L = Total length of each arc/link
			for (Pipe p : this.pipes.values()) {
				out.write("\n" + p.getStartNode().getNodeID() + "    " + p.getEndNode().getNodeID() + "    " + p.getLength());
			}
			out.write(";");

			out.write("\n\nparam E :=");  // NOTE: E = Elevation of each node
			for (Node n : this.nodes.values()) {
				// Selected based on network file analysis
				if (n.getNodeID() == this.source.getNodeID()) {
					// NOTE: It is necessary to use Head instead of Elevation when creating the network file
					//       so that the result of solver matches with the result of the old Jaltantra system
					//       for the input file "Sample_input_cycle_twoloop (Source Elevation changed) (acyclic).xls"
					//       which generates the network file that has md5sum "998a075a3545f6e8045a9c6538dbba2a".
					// NOTE: Constraint number 5 (i.e. "con5") of model "m1" and "m2" is the reason for this.
					//       It asks the solver to directly consider the head of the source to be equal to the
					//       elevation of the source. And, based on input file analysis, it was found that
					//       model files require the minimum pressure for source to be always 0. Hence, it is
					//       not possible to use true Elevation and "Minimum Pressure = Head - Elevation", and
					//       if we do so, then we get presolve error from the solvers.
					// NOTE: The formula "Pressure at each node = Head - Elevation" was found from:
					//       1. Note::getPressure()
					//       2. Analysis of:
					//          "Jaltantra website > Results section > Nodes tab > Elevation, Head and Pressure columns"
					out.write("\n" + n.getNodeID() + "   " + n.getHead());
					continue;
				}
				out.write("\n" + n.getNodeID() + "   " + n.getElevation());
			}
			out.write(";");

			out.write("\n\nparam P :=");  // NOTE: P = Minimum pressure required at each node
			for (Node n : this.nodes.values()) {
				// TODO: Check if we have to use getPressure() or getResidualPressure() method ?
				//       What is the difference between the two ?
				// out.write("\n" + n.getNodeID() + "   " + n.getPressure());
				out.write("\n" + n.getNodeID() + "   " + n.getResidualPressure());  // Selected based on network file analysis
			}
			out.write(";");

			out.write("\n\nparam D :=");  // NOTE: D = Demand of each node
			double totalDemand = 0.0;
			for (Node n : this.nodes.values()) {
				// NOTE: Demand of source node is always 0. Hence, we do not need to handle the source node separately.
				totalDemand += n.getDemand();
			}
			for (Node n : this.nodes.values()) {
				out.write("\n" + n.getNodeID() + "   " + n.getDemand());
			}
			out.write(";");

			out.write("\n\nparam d :=");  // NOTE: d = Diameter of each commercial pipe
			for (int i = 0; i < this.pipeCost.size(); ++i) {
				// NOTE: `this.pipeCost` is of type `ArrayList`, so `this.pipeCost.get(i)` will work in constant time
				out.write("\n" + i + "   " + this.pipeCost.get(i).getDiameter());
			}
			out.write(";");

			out.write("\n\nparam C :=");  // NOTE: C = Cost per unit length of each commercial pipe
			for (int i = 0; i < this.pipeCost.size(); ++i) {
				// NOTE: `this.pipeCost` is of type `ArrayList`, so `this.pipeCost.get(i)` will work in constant time
				out.write("\n" + i + "   " + this.pipeCost.get(i).getCost());
			}
			out.write(";");

			out.write("\n\nparam R :=");  // NOTE: R = Roughness of each commercial pipe
			for (int i = 0; i < this.pipeCost.size(); ++i) {
				// NOTE: `this.pipeCost` is of type `ArrayList`, so `this.pipeCost.get(i)` will work in constant time
				out.write("\n" + i + "   " + this.pipeCost.get(i).getRoughness());
			}
			out.write(";");

			out.write(String.format("\n\nparam Source := %d;\n", this.source.getNodeID()));  // NOTE: Source node ID
		} catch (IOException e) {
			loge("Network file creation failed due to IOException:");
			loge(e.getMessage());
			throw e;
		} finally {
			if (out != null) {
				out.close();
			}
		}
		// Network file successfully created

		// REFER: https://datacadamia.com/file/hash
		//        `Hash = Digest`
		// REFER: https://www.janbasktraining.com/community/sql-server/explain-the-differences-as-well-as-the-similarities-between-checksum-vs-hash
		//        According to me Hash and Checksum are also same in this context
		// REFER: https://www.geeksforgeeks.org/sha-256-hash-in-java/
		//        This is same as `shasum -a 256 FilePath`
		final String fileHash = findFileHashInHex(MessageDigest.getInstance("SHA-256"), new File(pathTo1FreshFile));

		final String hashedFileName = fileHash + ".R";
		final String pathTo2HashedFileName = SOLVER_ROOT_DIR + "/" + SOLVER_2_HASH_FILE_DIR + "/" + hashedFileName;
		if ((new File(pathTo2HashedFileName)).exists()) {
			logi("Request for this network was already submitted in the past: " +
				 "it may be in progress\uD83C\uDFC3 or finished\uD83C\uDFC1");
			// Delete this temporary file as this network was already given in the past
			if (!(new File(pathTo1FreshFile)).delete()) {
				loge("FIXME: Deletion of temporary network file failed: '" + pathTo1FreshFile + "'");
			}
			return "2-" + hashedFileName;
		}

		checkAndCreateDir(SOLVER_ROOT_DIR + "/" + SOLVER_2_HASH_FILE_DIR);

		// REFER: https://stackoverflow.com/questions/4645242/how-do-i-move-a-file-from-one-location-to-another-in-java
		boolean ok = (new File(pathTo1FreshFile)).renameTo(new File(pathTo2HashedFileName));
		if (!ok) {
			loge(String.format("FIXME: Failed to move the file: '%s' -> '%s'", pathTo1FreshFile, pathTo2HashedFileName));
			// Delete this temporary file because the rename operation failed
			if (!(new File(pathTo1FreshFile)).delete()) {
				loge("FIXME: Deletion of temporary network file failed: '" + pathTo1FreshFile + "'");
			}
			throw new Exception("Internal server error: problem with write permission");
		}

		logi(String.format("File renamed and moved successfully: '%s' -> '%s'", pathTo1FreshFile, pathTo2HashedFileName));

		logi("creating gams m1 and m2 model data file");
		String solverToGamsm1=SOLVER_ROOT_DIR + "/" + SOLVER_2_HASH_FILE_DIR + "/"+fileHash+"m1.gms";
		String solverToGamsm2=SOLVER_ROOT_DIR + "/" + SOLVER_2_HASH_FILE_DIR + "/"+fileHash+"m2.gms";
		createGamsModel(pathTo2HashedFileName,solverToGamsm1,"m1");
		createGamsModel(pathTo2HashedFileName,solverToGamsm2,"m2");
		logi("finished creating gams m1 and m2 model data file");
		return "0-" + hashedFileName;
	}

	/**
	 * Creates a GAMS model based on the provided input files and parameters.
	 *
	 * @param amplFilePath   The path to the AMPL input file.
	 * @param gamsFilePath   The path to the GAMS output file to be created.
	 * @param modelname      The name of the GAMS model.
	 * @throws IOException   If an I/O error occurs while reading the input file.
	 */
	public void createGamsModel(String amplFilePath,String gamsFilePath,String modelname) throws IOException {
		ArrayList<String> nodes = new ArrayList<>();
		ArrayList<String> pipes = new ArrayList<>();
		ArrayList<String> arcs = new ArrayList<>();
		ArrayList<String> arcsLength = new ArrayList<>();
		ArrayList<String> elevation = new ArrayList<>();
		ArrayList<String> pressure = new ArrayList<>();
		ArrayList<String> demand = new ArrayList<>();
		ArrayList<String> diameter = new ArrayList<>();
		ArrayList<String> pipeCost = new ArrayList<>();
		ArrayList<String> pipeRoughness = new ArrayList<>();
		String source="";

		File file = new File(amplFilePath);
		Scanner sc = new Scanner(file);

		int i = 1;
		while (sc.hasNext()) {
			if (i == 1) {
				readnodes(sc, nodes);
			} else if (i == 2) {
				readPipes(sc, pipes);
			} else if (i == 3) {
				readArcs(sc, arcs, arcsLength);
			} else if (i == 4) {
				readElevation(sc, elevation);
			} else if (i == 5) {
				readPressure(sc, pressure);
			} else if (i == 6) {
				readDemand(sc, demand);
			} else if (i == 7) {
				readDiameter(sc, diameter);
			} else if (i == 8) {
				readCost(sc, pipeCost);
			} else if (i == 9) {
				readRoughness(sc, pipeRoughness);
			} else {
				source=readHead(sc);
				break;
			}
			i++;
		}
		printSetsAndParametersAndModels(amplFilePath,gamsFilePath,nodes,pipes,arcs,arcsLength,diameter,elevation,pressure,demand,pipeCost,pipeRoughness,source,modelname);
		sc.close();
	}

	/**
	 * Prints the sets, parameters, and models to a GAMS output file based on the provided input and source.
	 *
	 * @param amplFilePath     The path to the AMPL input file.
	 * @param gamsFilePath     The path to the GAMS output file.
	 * @param nodes            The list of nodes.
	 * @param pipes            The list of pipes.
	 * @param arcs             The list of arcs.
	 * @param arcsLength       The list of arc lengths.
	 * @param diameter         The list of diameters.
	 * @param elevation        The list of elevations.
	 * @param pressure         The list of pressures.
	 * @param demand           The list of demands.
	 * @param pipeCost         The list of pipe costs.
	 * @param pipeRoughness    The list of pipe roughness values.
	 * @param source           The source information.
	 * @param modelName        The name of the model.
	 * @throws IOException     If an I/O error occurs while writing to the output file.
	 */
	void printSetsAndParametersAndModels(String amplFilePath, String gamsFilePath, ArrayList<String> nodes, ArrayList<String> pipes, ArrayList<String> arcs, ArrayList<String> arcsLength, ArrayList<String> diameter, ArrayList<String> elevation, ArrayList<String> pressure, ArrayList<String> demand, ArrayList<String> pipeCost, ArrayList<String> pipeRoughness, String source, String modelName) throws IOException {

		BufferedWriter out = null;

		try {
			FileWriter fstream = new FileWriter(gamsFilePath, true); //true tells to append data.
			out = new BufferedWriter(fstream);

			printSets(out,nodes,pipes,arcs,source);
			printParameters(out,arcsLength,diameter,elevation,pressure,demand,pipeCost,pipeRoughness,source);
			if(modelName.equals("m1")){
				printModelm1(out);
			}
			else printModelm2(out);
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

	/**
	 * Writes the model M2 equations and solve statement to the specified output writer.
	 *
	 * @param out The output writer to write the model M2 equations.
	 * @throws IOException If an I/O error occurs while writing to the output writer.
	 */
	void printModelm2(BufferedWriter out) throws IOException {
		out.write("Scalar omega  /10.68/;\n");
		out.write("Scalar bnd ;\n");
		out.write("Scalar qm;\n");
		out.write("Scalar q_M;\n");

		out.write("\n");

		out.write("bnd = sum(src,D(src));\n");
		out.write("q_M=-bnd;\n");
		out.write("qm=0;\n");

		out.write("\n");

		out.write("Variable l(nodes,j,pipes); \n");
		out.write("l.lo(nodes,j,pipes)= 0;\n");

		out.write("\n");

		out.write("Variable q1(nodes,j);\n");
		out.write("q1.lo(nodes,j)=qm;\n");
		out.write("q1.up(nodes,j)=q_M;\n");

		out.write("\n");

		out.write("Variable q2(nodes,j);\n");
		out.write("q2.lo(nodes,j)=qm;\n");
		out.write("q2.up(nodes,j)=q_M;\n");

		out.write("\n");

		out.write("Variables z;\n");

		out.write("\n");

		out.write("Variable h(nodes);\n");

		out.write("\n");

		out.write("Equations cost \"objective function\",bound1(nodes,j,pipes),cons1(nodes),cons2(nodes),cons3(nodes,j),cons5(src), cons4(nodes,j), cons6(nodes,j) ;\n");

		out.write("cost..  z=e=sum(arcs(nodes,j),sum(pipes,l(arcs,pipes)*c(pipes)));\n");

		out.write("bound1(nodes,j,pipes)$arcs(nodes,j).. l(nodes,j,pipes) =l= Len(nodes,j);\n");
		out.write("cons1(nodes).. sum(arcs(j,nodes),(q1(arcs)-q2(arcs))) =e= sum(arcs(nodes,j),(q1(arcs)-q2(arcs))) + D(nodes);\n");
		out.write("cons2(nodes).. h(nodes) =g= E(nodes) + P(nodes);\n");
		out.write("cons3(arcs(nodes,j)).. h(nodes)-h(j)=e=sum(pipes,(((q1(arcs)*0.001)**1.852 - (q2(arcs)*0.001)**1.852)*omega*l(arcs,pipes))/((R(pipes)**1.852)*(dis(pipes)/1000)**4.87));\n");
		out.write("cons4(arcs(nodes,j)).. sum(pipes,l(arcs,pipes)) =e=Len(arcs);\n");
		out.write("cons5(src)..  h(src)=e= sum(srcs,E(srcs));\n");
		out.write("cons6(arcs(nodes,j)).. q1(arcs)*q2(arcs) =l= q_M*qm;\n");

		out.write("\n");

		out.write("model m2  /all/  ;\n");
		// System.out.println("Option threads=4;\n");
		// System.out.println("m2.optfile =1;\n");
		out.write("solve m2 using minlp minimizing z ;\n");
	}

	/**
	 * Writes the model M1 equations and solve statement to the specified output writer.
	 *
	 * @param out The output writer to write the model M1 equations.
	 * @throws IOException If an I/O error occurs while writing to the output writer.
	 */
	private void printModelm1(BufferedWriter out) throws IOException {

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
// 		out.write("m1.optfile =1;\n");
		out.write("solve m1 using minlp minimizing z ;\n");

	}

	/**
	 * Writes the parameters to the specified output writer.
	 *
	 * @param out             The output writer to write the parameters.
	 * @param arcsLength      The list of arcs lengths.
	 * @param diameter        The list of pipe diameters.
	 * @param elevation       The list of node elevations.
	 * @param pressure        The list of node pressures.
	 * @param demand          The list of node demands.
	 * @param pipeCost        The list of pipe costs.
	 * @param pipeRoughness   The list of pipe roughness values.
	 * @param source          The source node.
	 * @throws IOException   If an I/O error occurs while writing to the output writer.
	 */
	private void printParameters(BufferedWriter out, ArrayList<String> arcsLength, ArrayList<String> diameter, ArrayList<String> elevation, ArrayList<String> pressure, ArrayList<String> demand, ArrayList<String> pipeCost, ArrayList<String> pipeRoughness, String source) throws IOException {

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

	/**
	 * Writes the sets to the specified output writer.
	 *
	 * @param out        The output writer to write the sets.
	 * @param nodes      The list of nodes.
	 * @param pipes      The list of pipes.
	 * @param arcs       The list of arcs.
	 * @param source     The source node.
	 * @throws IOException   If an I/O error occurs while writing to the output writer.
	 */
	private void printSets(BufferedWriter out, ArrayList<String> nodes, ArrayList<String> pipes, ArrayList<String> arcs, String source) throws IOException {

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


	// printSetsAndParametersAndModels();
	/**
	 * Reads nodes from the specified scanner and adds them to the given list of nodes.
	 *
	 * @param sc     The scanner to read input from.
	 * @param nodes  The list to store the nodes.
	 */
	public void readnodes(Scanner sc,ArrayList<String> nodes){
		String inp=sc.nextLine();
		String str[]=inp.split("\\s+");
		String node="";
		for(int i=3;i<str.length-1;i++)
			node=node+str[i]+", ";
		node=node+str[str.length-1].substring(0,str[str.length-1].length()-1);
		nodes.add(node);
		sc.nextLine();
//		    System.out.println(nodes);
	}

	/**
	 * Reads pipes from the specified scanner and adds them to the given list of pipes.
	 *
	 * @param sc     The scanner to read input from.
	 * @param pipes  The list to store the pipes.
	 */
	public void readPipes(Scanner sc,ArrayList<String> pipes){
		String inp=sc.nextLine();
		String str[]=inp.split("\\s+");
		String pipe="";
		for(int i=3;i<str.length-1;i++)
			pipe=pipe+str[i]+", ";
		pipe=pipe+str[str.length-1].substring(0,str[str.length-1].length()-1);
		pipes.add(pipe);
		sc.nextLine();
//		    System.out.println(pipes);
	}

	/**
	 * Reads the arcs and their lengths from the specified scanner and populates the given lists.
	 *
	 * @param sc          The scanner to read input from.
	 * @param arcs        The list to store the arcs.
	 * @param arcsLength  The list to store the lengths of arcs.
	 */
	public void readArcs(Scanner sc,ArrayList<String> arcs,ArrayList<String> arcsLength){
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
		String len=str[0]+"    ."+str[1]+"    "+str[2].substring(0,str[2].length()-1);
		arcsLength.add(len);
//		    System.out.println(pipes);
		sc.nextLine();
	}

	/**
	 * Reads diameter from the specified scanner and adds them to the given list of diameter.
	 *
	 * @param sc     The scanner to read input from.
	 * @param diameter  The list to store the diameter.
	 */
	public void readDiameter(Scanner sc, ArrayList<String> diameter){
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
		String dia=str[0]+" "+str[1].substring(0,str[1].length()-1);
		diameter.add(dia);
		sc.nextLine();
	}

	/**
	 * Reads elevation from the specified scanner and adds them to the given list of elevation.
	 *
	 * @param sc     The scanner to read input from.
	 * @param elevation  The list to store the elevation.
	 */
	public void readElevation(Scanner sc, ArrayList<String> elevation){
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
		String ele=str[0]+" "+str[1].substring(0,str[1].length()-1);
		elevation.add(ele);
		sc.nextLine();
	}

	/**
	 * Reads pressure from the specified scanner and adds them to the given list of pressure.
	 *
	 * @param sc     The scanner to read input from.
	 * @param pressure  The list to store the pressure.
	 */
	public void readPressure(Scanner sc, ArrayList<String> pressure){
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
		String pres=str[0]+" "+str[1].substring(0,str[1].length()-1);
		pressure.add(pres);
		sc.nextLine();
	}

	/**
	 * Reads demand from the specified scanner and adds them to the given list of demand.
	 *
	 * @param sc     The scanner to read input from.
	 * @param demand  The list to store the demand.
	 */
	public void readDemand(Scanner sc, ArrayList<String> demand){
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
		String dem=str[0]+" "+str[1].substring(0,str[1].length()-1);
		demand.add(dem);
		sc.nextLine();
	}

	/**
	 * Reads pipeCost from the specified scanner and adds them to the given list of pipeCost.
	 *
	 * @param sc     The scanner to read input from.
	 * @param pipeCost  The list to store the pipeCost.
	 */
	public void readCost(Scanner sc, ArrayList<String> pipeCost){
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
		String cost=str[0] + " " + str[1].substring(0,str[1].length()-1);
		pipeCost.add(cost);
		sc.nextLine();
	}

	/**
	 * Reads pipeRoughness from the specified scanner and adds them to the given list of pipeRoughness.
	 *
	 * @param sc     The scanner to read input from.
	 * @param pipeRoughness  The list to store the pipeRoughness.
	 */
	public void readRoughness(Scanner sc, ArrayList<String> pipeRoughness){
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
		String roughness=str[0]+" "+str[1].substring(0,str[1].length()-1);
		pipeRoughness.add(roughness);
		sc.nextLine();
	}

	/**
	 * Reads the header from the specified scanner and returns the extracted source.
	 *
	 * @param sc The scanner to read input from.
	 * @return The extracted source.
	 */
	public String readHead(Scanner sc){
		String inp=sc.nextLine();
		String source=inp.split("\\s+")[3];
		source=source.substring(0,source.length()-1);
		return source;
	}


	/**
	 * Launch "CalculateNetworkCost_JaltantraLauncher.sh" and wait for 2 seconds. If the process
	 * exits within 2 seconds, then log it (saying that we need to check into this matter) and return.
	 *
	 * @param networkFilePath path to the network file which is the be solved
	 * @throws Exception if `SOLVER_ROOT_DIR` directory does not exist,
	 *                   or function is unable to start the launcher process,
	 *                   or exception is thrown by `process.wairFor(...)`
	 */
	private void launchCalculateNetworkCost(final String networkFilePath) throws Exception {
		if (!(new File(SOLVER_ROOT_DIR)).exists()) {
			loge("FIXME: SOLVER_ROOT_DIR directory does not exist: ' " + SOLVER_ROOT_DIR + "'");
			throw new Exception("Internal server error: SOLVER_ROOT_DIR directory does not exist");
		}

		// REFER: https://mkyong.com/java/how-to-execute-shell-command-from-java/
		// REFER: https://www.geeksforgeeks.org/how-to-execute-native-shell-commands-from-java-program/
		try {
			logd(System.getProperty("os.name"));

			// NOTE: This process cannot be run on Windows
			ProcessBuilder pb = new ProcessBuilder(
					"bash",
					SOLVER_ROOT_DIR + "/CalculateNetworkCost_JaltantraLauncher.sh",
					networkFilePath,
					SOLVER_EXECUTION_TIME,
					run_Time
			);
			pb.directory(new File(System.getProperty("user.home")));
			Process process = pb.start();

			boolean exited = process.waitFor(2, TimeUnit.SECONDS);
			if (exited)
				logi("CHECKME: CalculateNetworkCost.py exited within 2 seconds, probably due to some issue");
		} catch (IOException | InterruptedException e) {
			e.printStackTrace();
		}
	}

	/**
	 * Check the execution status of "CalculateNetworkCost.py"
	 *
	 * @return int
	 * <li>-1 if "pathToNetworkSpecificDirectory/0_status" file does not exist</li>
	 * <li>-2 if solver status has unknown/unhandled value (refer code and
	 * documentation of "CalculateNetworkCost.py", and update this method)</li>
	 * <li>1 if solver is "running"</li>
	 * <li>2 if solver has "finished"</li>
	 * <li>3 if solver has "successfully finished"</li>
	 * @throws Exception if "launch_error" occurred
	 */
	private int checkSolverResultStatus(final String pathToNetworkSpecificDirectory) throws Exception {
		// REFER: CalculateNetworkCost.py docs for this
		if (!(new File(pathToNetworkSpecificDirectory + "/0_status")).exists())
			return -1;
		// NOTE: We use strip() as a safety measure to remove useless spaces
		// NOTE: We have to use startsWith() method because of trailing new line character(s) and other extra
		//       things which "CalculateNetworkCost.py" may be writing to the "0_status" file
		final String status = readFileAsString(pathToNetworkSpecificDirectory + "/0_status").strip();
		if (status.startsWith("launch_error")) {
			loge("FIXME: Failed to start the solver for network file: " + pathToNetworkSpecificDirectory.substring(pathToNetworkSpecificDirectory.lastIndexOf("/") + 1));
			throw new Exception("Internal server error: failed to start the solver for this network =" + pathToNetworkSpecificDirectory.substring(pathToNetworkSpecificDirectory.lastIndexOf("/") + 1) + "\n" + status.substring(status.indexOf("\n") + 1));
		}
		if (status.startsWith("running"))
			return 1;
		if (status.startsWith("finished"))
			return 2;
		if (status.startsWith("success"))
			return 3;
		return -2;
	}

	/**
	 * Extract solution for a given network and store it in:
	 * <br>&emsp;• resultPipes
	 * <br>&emsp;• resultCost
	 * <br>&emsp;• resultPumps
	 *
	 * @param pathToNetworkSpecificDirectory path to the directory which was created by CalculateNetworkCost.py
	 *                                       to store all the data related to a network file given to it to solve
	 * @return true if solver result was successfully extracted
	 * @throws Exception if any one of the below condition occurs:
	 *                   <br>&emsp;• result file does not exist
	 *                   <br>&emsp;• failure to read the result file
	 *                   <br>&emsp;• solver failed to find the solution due to some reason
	 */
	private boolean extractSolverResult(final String pathToNetworkSpecificDirectory) throws Exception {
		final String pathToResultFile = pathToNetworkSpecificDirectory + "/0_result.txt";
		if (!(new File(pathToResultFile)).exists())
			throw new Exception("Internal server error: result file does not exist");

		final String result = readFileAsString(pathToResultFile);
		final String[] resultLines = result.split("\n");
		logi("0_result.txt first line = " + resultLines[0]);
		if (resultLines[0].equals("False")) {
			loge("CHECKME: Solvers failed to find the solution for the network with hash = " + pathToNetworkSpecificDirectory.substring(pathToNetworkSpecificDirectory.lastIndexOf("/") + 1));
			throw new Exception("Internal server error: solvers failed to find the solution");
		}
		logi("Best result found for network with hash = " + pathToNetworkSpecificDirectory.substring(pathToNetworkSpecificDirectory.lastIndexOf("/") + 1));
		logi("Solver = " + resultLines[1]);
		logi("Model = " + resultLines[2]);
		logi("Best result = " + resultLines[4]);

		// REFER: https://www.geeksforgeeks.org/how-to-execute-native-shell-commands-from-java-program/
		ArrayList<String> extractedSolutionLines = new ArrayList<>();
		TreeMap<PipeCost, Double> cumulativePipeLength = new TreeMap<>();
		this.resultPipes = new ArrayList<>();
		this.resultCost = new ArrayList<>();
		this.resultPumps = new ArrayList<>();

		try {
			// This process cannot be run on Windows
			ProcessBuilder pb = new ProcessBuilder(
					"python3",
					SOLVER_ROOT_DIR + "/CalculateNetworkCost_ExtractResultFromAmplOutput.py",
					resultLines[3],
					pathToNetworkSpecificDirectory + "/0_graph_network_data_testcase.R",
					"1"
			);
			pb.directory(new File(System.getProperty("user.home")));
			Process process = pb.start();

			BufferedReader reader = new BufferedReader(new InputStreamReader(process.getInputStream()));

			logd("*** The solution is ***");
			String line;
			while ((line = reader.readLine()) != null) {
				logd(line);
				extractedSolutionLines.add(line);
			}

			int exitVal = process.waitFor();
			if (exitVal != 0) {
				logi("extractSolverResult(...): exitVal = 0");
				loge("********************");
				loge("Log of CalculateNetworkCost_ExtractResultFromAmplOutput.py");
				loge("********************");
				BufferedReader readerErr = new BufferedReader(new InputStreamReader(process.getErrorStream()));
				while ((line = readerErr.readLine()) != null) {
					loge(line);
				}
				loge("********************");
				return false;
			}

			int solutionLineIdx = 0;
			int numberOfHeads = Integer.parseInt(extractedSolutionLines.get(solutionLineIdx));
			solutionLineIdx += 1;
			for (int iHead = 0; iHead < numberOfHeads; ++iHead) {
				String[] arr = extractedSolutionLines.get(solutionLineIdx).split(" ");
				solutionLineIdx += 1;

				int node = Integer.parseInt(arr[0]);
				double head = Double.parseDouble(arr[1]);

				this.nodes.get(node).setHead(head);
			}

			int numberOfFlows = Integer.parseInt(extractedSolutionLines.get(solutionLineIdx));
			solutionLineIdx += 1;
			for (int iFlow = 0; iFlow < numberOfFlows; ++iFlow) {
				String[] arr = extractedSolutionLines.get(solutionLineIdx).split(" ");
				solutionLineIdx += 1;

				int arcSourceVertex = Integer.parseInt(arr[0]);
				int arcDestinationVertex = Integer.parseInt(arr[1]);
				double flow = Double.parseDouble(arr[2]);
				if (flow < 0.0) {
					logi("IMPORTANT: -ve flow indicates that water flows from destination to source vertex");
					// logi("IMPORTANT: Swapping source and destination vertex as flow is -ve");
					// int temp = arcDestinationVertex;
					// arcDestinationVertex = arcSourceVertex;
					// arcSourceVertex = temp;
				}

				Pipe lCurrArc = null;
				for (Pipe p : this.getPipes().values()) {
					if (p.getStartNode().getNodeID() != arcSourceVertex || p.getEndNode().getNodeID() != arcDestinationVertex)
						continue;
					lCurrArc = p;
					break;
				}
				if (lCurrArc == null) {
					loge("arcSourceVertex = " + arcSourceVertex + ", arcDestinationVertex = " + arcDestinationVertex);
					throw new Exception("Internal server error: Failed to extract the solution");
				}
				lCurrArc.setFlow(flow);
			}

			int numberOfArcs = Integer.parseInt(extractedSolutionLines.get(solutionLineIdx));
			solutionLineIdx += 1;
			for (int arcNum = 0; arcNum < numberOfArcs; ++arcNum) {
				String[] arr = extractedSolutionLines.get(solutionLineIdx).split(" ");
				solutionLineIdx += 1;
				int arcSourceVertex = Integer.parseInt(arr[0]);
				int arcDestinationVertex = Integer.parseInt(arr[1]);
				int arcOptimalPipeCount = Integer.parseInt(arr[2]);
				Pipe lCurrArc = null;
				for (Pipe p : this.getPipes().values()) {
					if (p.getStartNode().getNodeID() != arcSourceVertex || p.getEndNode().getNodeID() != arcDestinationVertex)
						continue;
					lCurrArc = p;
					break;
				}
				if (lCurrArc == null) {
					throw new Exception("Internal server error: Failed to extract the solution");
				}
				for (int pipeNum = 0; pipeNum < arcOptimalPipeCount; ++pipeNum) {
					// Prefix "l" denotes local variable
					String[] lArr2 = extractedSolutionLines.get(solutionLineIdx).split(" ");
					solutionLineIdx += 1;
					int lPipeIdx = Integer.parseInt(lArr2[0]);
					double lPipeLength = Double.parseDouble(lArr2[1]);

					PipeCost lCurrPipeCost = this.getPipeCost().get(lPipeIdx);
					double lPipeDiameter = lCurrPipeCost.getDiameter();
					double lPipeFlow = lCurrArc.getFlow();
					double lPipeHeadLoss = Util.HWheadLoss(lPipeLength, lPipeFlow, lCurrArc.getRoughness(), lPipeDiameter);
					boolean lPipePressureExceeded = false;
					if (generalProperties.max_pipe_pressure > 0.0)
						if (lCurrArc.getStartNode().getPressure() > generalProperties.max_pipe_pressure || lCurrArc.getEndNode().getPressure() > generalProperties.max_pipe_pressure)
							lPipePressureExceeded = true;

					int linkId = -1;
					for (Pipe p : this.pipes.values()) {
						if (arcSourceVertex == p.getStartNode().getNodeID()
							&& arcDestinationVertex == p.getEndNode().getNodeID()) {
							linkId = p.getPipeID();
							break;
						}
					}
					if (linkId == -1) {
						loge(String.format(
								"Link ID could not be found for " +
								"arcSourceVertex=%d, arcDestinationVertex=%d, lPipeLength=%f, lPipeDiameter=%f",
								arcSourceVertex,
								arcDestinationVertex,
								lPipeLength,
								lPipeDiameter
						));
					}
					this.resultPipes.add(new PipeStruct(
							linkId, arcSourceVertex, arcDestinationVertex, lPipeLength, lPipeDiameter,
							lCurrPipeCost.getRoughness(), lPipeFlow, lPipeHeadLoss, lPipeHeadLoss * 1000 / lPipeLength,
							Util.waterSpeed(lPipeFlow, lPipeDiameter), lPipeLength * lCurrPipeCost.getCost(),
							false, lPipePressureExceeded, lCurrArc.getFlowchoice() == FlowType.PRIMARY,
							lCurrArc.getPumpHead(), lCurrArc.getPumpPower(), lCurrArc.getValveSetting()
					));
					if (pumpGeneralProperties.pump_enabled && lCurrArc.getPumpPower() > 0) {
						double presentvaluefactor = Util.presentValueFactor(pumpGeneralProperties.discount_rate, pumpGeneralProperties.inflation_rate, pumpGeneralProperties.design_lifetime);
						double primarycoeffecient = 365 * presentvaluefactor * generalProperties.supply_hours * pumpGeneralProperties.energycost_per_kwh * pumpGeneralProperties.energycost_factor;
						double secondarycoeffecient = esrGeneralProperties.esr_enabled ? 365 * presentvaluefactor * esrGeneralProperties.secondary_supply_hours * pumpGeneralProperties.energycost_per_kwh * pumpGeneralProperties.energycost_factor : 0;
						double power = lCurrArc.getPumpPower();

						double energycost = power * (lCurrArc.getFlowchoice() == FlowType.PRIMARY ? primarycoeffecient : secondarycoeffecient);
						double capitalcost = power * pumpGeneralProperties.capitalcost_per_kw;


						resultPumps.add(new ResultPumpStruct(lCurrArc.getPipeID(),
								lCurrArc.getPumpHead(),
								power,
								energycost,
								capitalcost,
								energycost + capitalcost));
					}
					if (!cumulativePipeLength.containsKey(lCurrPipeCost))
						cumulativePipeLength.put(lCurrPipeCost, lPipeLength);
					else
						cumulativePipeLength.put(lCurrPipeCost, cumulativePipeLength.get(lCurrPipeCost) + lPipeLength);
				}
			}

			double cumulativeCost = 0;
			for (Map.Entry<PipeCost, Double> entry : cumulativePipeLength.entrySet()) {
				double cost = 0;
				cost = entry.getValue() * entry.getKey().getCost();
				cumulativeCost += cost;
				CommercialPipeStruct resultcommercialPipe = new CommercialPipeStruct(
						entry.getKey().getDiameter(),
						cost,
						entry.getValue(),
						cumulativeCost,
						entry.getKey().getRoughness()
				);
				resultCost.add(resultcommercialPipe);
			}
			return true;
		} catch (IOException | InterruptedException e) {
			loge(e.toString());
			e.printStackTrace();
		}

		return false;
	}

	// ---

	public HashMap<Integer, Node> getNodes() {
		return nodes;
	}

	public HashMap<Integer, Pipe> getPipes() {
		return pipes;
	}

	public List<PipeCost> getPipeCost() {
		return pipeCost;
	}

}