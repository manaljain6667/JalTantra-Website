package optimizer;

import com.google.ortools.linearsolver.MPSolver.ResultStatus;
import optimizer.Pipe.FlowType;
import structs.*;

import java.io.*;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.security.MessageDigest;
import java.text.SimpleDateFormat;
import java.text.*;
import java.util.*;
import java.util.concurrent.TimeUnit;


// This class is responsible for optimization of the network
public class Optimizer {
	public static void logd(String str) {
		System.err.println("FM: DEBUG: " + str);
	}

	public static void logi(String str) {
		System.err.println("FM: INFO: " + str);
	}

	public static void logw(String str) {
		System.err.println("FM: WARNING: " + str);
	}

	public static void loge(String str) {
		System.err.println("FM: ERROR: " + str);
	}

	public ArrayList<PipeStruct> resultPipes;
	public ArrayList<CommercialPipeStruct> resultCost;
	public ArrayList<ResultPumpStruct> resultPumps;

	private HashMap<Integer,Node> nodes; // nodeID,node
	private HashMap<Integer,Pipe> pipes; // pipeid,pipe
	private List<PipeCost> pipeCost;
	private List<EsrCost> esrCost;
	
	//following three are data structures received from server request
	private GeneralStruct generalProperties;
	private EsrGeneralStruct esrGeneralProperties;
	private PumpGeneralStruct pumpGeneralProperties;
	
	private Node source;	//source node of the network
	private double totalDemand;	//total demand of the network in litres/sec
	
	private Problem problem;	//contains the ILP model
	
	//helper strings for exporting the ILP model 
	private StringBuilder lpmodelstring = new StringBuilder();
	private StringBuilder lpmodelvar = new StringBuilder();
	
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
	private double minEsrHeight = 0;
	private double maxEsrHeight = 25;
	
	// minimum and maximum power allowed for a pump
	private double minPumpPower = 1;
	private double maxPumpPower = 10000;
	
	// maximum pressure head that can be provided by a pump
	private int maxPumpHead = 10000;
	
	//string used for EPANET output file
	private String coordinatesString;
	
	//container for pump and valve information
	private PumpManualStruct[] pumpManualArray;
	private ValveStruct[] valves;
	
	//create an instance of optimizer  
	public Optimizer(NodeStruct[] nodeStructs, PipeStruct[] pipeStructs, CommercialPipeStruct[] commercialPipeStructs, GeneralStruct generalStruct, EsrGeneralStruct esrGeneralProperties, EsrCostStruct[] esrCostsArray, PumpGeneralStruct pumpGeneralProperties, PumpManualStruct[] pumpManualArray, ValveStruct[] valves) throws Exception{
		nodes = new HashMap<Integer, Node>();
		pipes = new HashMap<Integer, Pipe>();
		pipeCost = new ArrayList<PipeCost>();
		esrCost = new ArrayList<EsrCost>();
		coordinatesString = "";
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
		for(NodeStruct node : nodeStructs){
			double minPressure = node.minpressure == 0 ? generalProperties.min_node_pressure : node.minpressure;
			Node n = new Node(node.elevation, node.demand, node.nodeid, minPressure, node.nodename, 24/generalProperties.supply_hours, usedNodeIDs);
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
		for(PipeStruct pipe : pipeStructs){
			double roughness = pipe.roughness == 0 ? generalProperties.def_pipe_roughness : pipe.roughness;
			
			Node startNode = nodes.get(pipe.startnode);
            if(startNode==null){
            	throw new Exception("Invalid startNode:" + pipe.startnode + " provided for pipe ID:"+pipe.pipeid);
            }
            
            Node endNode = nodes.get(pipe.endnode);
            if(endNode==null){
            	throw new Exception("Invalid endNode:" + pipe.endnode + " provided for pipe ID:"+pipe.pipeid);
            }
			
			Pipe p = new Pipe(pipe.length, startNode, endNode, pipe.diameter, roughness, pipe.pipeid, pipe.parallelallowed, usedPipeIDs);
			usedPipeIDs.add(p.getPipeID());
			pipes.put(p.getPipeID(), p);
		}
		
		//initialize the commercial pipe information
		for(CommercialPipeStruct commercialPipe : commercialPipeStructs){
			double roughness = commercialPipe.roughness == 0 ? generalProperties.def_pipe_roughness : commercialPipe.roughness;
			pipeCost.add(new PipeCost(commercialPipe.diameter, commercialPipe.cost, Double.MAX_VALUE, roughness));
		}
		this.valves = valves;
		
		//default model number is 0 for only pipe optimization
		modelNumber = 0;
		
		//if ESR optimization enabled, initialize ESR properties and set modelnumber
		if(esrGeneralProperties!=null && esrGeneralProperties.esr_enabled){
			this.esrGeneralProperties = esrGeneralProperties;
			
			if(esrGeneralProperties.secondary_supply_hours==0){
				throw new Exception("ESR option is enabled, but secondary supply hours is provided as zero.");
			}
			
			if(esrGeneralProperties.esr_capacity_factor==0){
				throw new Exception("ESR option is enabled, but esr capacity factor is provided as zero.");
			}
			
			secondaryFlowFactor = generalProperties.supply_hours/esrGeneralProperties.secondary_supply_hours;
			esrCapacityFactor = esrGeneralProperties.esr_capacity_factor;
			maxEsrHeight = esrGeneralProperties.max_esr_height;
			
			modelNumber = 9;
			
			for(EsrCostStruct esrcost : esrCostsArray){
				esrCost.add(new EsrCost(esrcost.mincapacity,
										esrcost.maxcapacity,
										esrcost.basecost,
										esrcost.unitcost));
			}
		}
		
		//if pump enabled, initialize pump properties
		if(pumpGeneralProperties!=null && pumpGeneralProperties.pump_enabled){
			this.pumpGeneralProperties = pumpGeneralProperties;
			this.pumpManualArray = pumpManualArray;
			this.minPumpPower = pumpGeneralProperties.minpumpsize;
			
			if(pumpGeneralProperties.efficiency==0)
				throw new Exception("Pump option is enabled, but pump efficiency is provided as zero.");
			
			if(pumpGeneralProperties.design_lifetime==0)
				throw new Exception("Pump option is enabled, but design lifetime is provided as zero.");
		}
		
		//set total demand required for the network
		totalDemand = getTotalCapacity();		
	}

	/**
	 * Check whether the network structure is valid or not
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

	//get the total water supply that flows through a node in litres per second
	//recursively compute for all downstream nodes
	//simultaneously set the flow through outgoing pipes 
	//simultaneously set the downstreamnodes property (only demand nodes considered)
	private double getNodeSupply(Node node){
		double sum=node.getDemand();
		double supply;
		for(Pipe pipe : node.getOutgoingPipes()){
			Node e = pipe.getEndNode();
			supply = getNodeSupply(e);
			pipe.setFlow(supply);
			sum += supply;
			
			node.addToDownstreamNodes(e.getDownstreamNodes());
			if(e.getDemand()!=0)
				node.addToDownstreamNodes(e);
		}
		return sum;
	}
	
	//return the total ESR capacity required in the network in litres 
	private double getTotalCapacity(){
		double sum = 0;
		for(Node n: nodes.values()){
			sum = sum + n.getRequiredCapacity(esrCapacityFactor);
		}
		return sum;
	}
	
	//get the total water supply that flows through a node in litres per second
	//recursively compute for all downstream nodes
	//simultaneously set the flow through outgoing pipes 
	//simultaneously set the downstreamnodes property (includes zero demand nodes as well)	
	private static double getNodeSupply_gen(Node node){
		double sum=node.getDemand();
		double supply;
		for(Pipe pipe : node.getOutgoingPipes()){
			Node e = pipe.getEndNode();
			supply = getNodeSupply_gen(e);
			pipe.setFlow(supply);
			sum += supply;
			
			node.addToDownstreamNodes(e.getDownstreamNodes());
			node.addToDownstreamNodes(e);
		}
		return sum;
	}
	
	//sets the sourcetonodepipes property of a node
	//recursively call this for all downstream nodes
	//simultaneously set the upstreamnodes property of the node (only considering nodes with demand)
	private static void setSourceToNodePipes(Node node){
		Node n;
		for(Pipe pipe : node.getOutgoingPipes()){
			n = pipe.getEndNode();
			n.addToSourceToNodePipes(node.getSourceToNodePipes());
			n.addToSourceToNodePipes(pipe);
			
			n.addToUpstreamNodes(node.getUpstreamNodes());
			if(node.getDemand()!=0)
				n.addToUpstreamNodes(node);
			
			setSourceToNodePipes(n);
		}
	}
	
	//sets the sourcetonodepipes property of a node
	//recursively call this for all downstream nodes
	//simultaneously set the upstreamnodes property of the node (includes node with zero demand)		
	private static void setSourceToNodePipes_gen(Node node){
		Node n;
		for(Pipe pipe : node.getOutgoingPipes()){
			n = pipe.getEndNode();
			n.addToSourceToNodePipes(node.getSourceToNodePipes());
			n.addToSourceToNodePipes(pipe);
			
			n.addToUpstreamNodes(node.getUpstreamNodes());
			n.addToUpstreamNodes(node);
			
			setSourceToNodePipes_gen(n);
		}
	}
	
	//set the objective cost of the ILP 
	//includes capital cost of pipes and capital+energy cost of pumps(if enbaled)
	private void setObjectiveCost() throws Exception{
		Linear linear = new Linear();
		int j=0;
		for(Pipe pipe : pipes.values()){
			int i = pipe.getPipeID();		
			j=0;
			for(PipeCost entry : pipeCost){	
				problem.setVarType("l_"+i+"_"+j, Double.class);
				problem.setVarLowerBound("l_"+i+"_"+j, 0);
				if(pipe.getDiameter()==0) // cost contributes only if diameter is to be computed
					linear.add(entry.getCost(), "l_"+i+"_"+j);
				if(pipe.isAllowParallel()){
					problem.setVarType("p_"+i+"_"+j, Boolean.class);
					linear.add(entry.getCost() * pipe.getLength(), "p_"+i+"_"+j);
				}
				j++;
			}
			if(pipe.isAllowParallel())
				problem.setVarType("p_"+i, Boolean.class);	
			
			if(pumpGeneralProperties.pump_enabled){
				problem.setVarType("pumphead_"+i, Double.class);
				problem.setVarLowerBound("pumphead_"+i, 0);
				problem.setVarUpperBound("pumphead_"+i, maxPumpHead);
							
				problem.setVarType("pumppower_"+i, Double.class);
				problem.setVarLowerBound("pumppower_"+i, 0);
								
				problem.setVarType("pumphelper_"+i, Boolean.class);
				
				double presentvaluefactor = Util.presentValueFactor(pumpGeneralProperties.discount_rate, pumpGeneralProperties.inflation_rate, pumpGeneralProperties.design_lifetime);
				double primarycoefficient = presentvaluefactor*365*generalProperties.supply_hours*pumpGeneralProperties.energycost_per_kwh;
								
				//capital + energy cost
				linear.add(pumpGeneralProperties.capitalcost_per_kw + primarycoefficient, "pumppower_"+i);
			}
			
		}
		
		//problem.setObjective(linear, OptType.MIN);
		problem.setObjective(linear, false);
	}
	
	//add variables to the ILP
	private void addVariables() throws Exception{
		// for each link i, for each commercial pipe j:
		// l_i_j : length of commercial pipe j of link i
		// p_i_j : boolean if ith link has jth pipe diameter in parallel
		// p_i : if ith link has no parallel pipe despite being allowed to have one
		// f_i denotes whether flow in pipe is primary = 1 or secondary = 0
		
		// helper variables to linearize product of binary and continous functions:
		// y_i_j = f_i * l_i_j
		// yp_i_j = p_i_j * f_i
		
		// head_i is the incoming head at node i
		// introduce head_i_j : i is the start node, j is the id of the pipe
		// head_i_j is the head provided by i to the outgoing pipe j
		// if there is no ESR at i or if there is an ESR at i but water to pipe j still comes from the primary network, head_i_j = head_i
		// if there is an ESR at i and it is responsible for providing water to the pipe j, head_i_j = elevation_i + esr_i 
		// i.e. head provided is determined by the esr height and the location elevation
		// k_i_j = s_i_i && !f_j
		// head_i_j = k_i_j * (elevation_i + esr_i - head_i) + head_i
		// ehead_i = elevation_i + esr_i - head_i
		// yhead_i_j = k_i_j * ehead_i
		for(Pipe pipe : pipes.values()){
			int i = pipe.getPipeID();	
			problem.setVarType("f_"+i, Boolean.class);
			
			problem.setVarType("headloss_"+i, Double.class);
			problem.setVarLowerBound("headloss_"+i, 0);
			
			if(pipe.isAllowParallel()){
				problem.setVarType("p_"+i, Boolean.class);	
				problem.setVarType("yp_"+i, Boolean.class);
			}
			
			Node startNode = pipe.getStartNode();
			if(startNode.getDemand()>0){
				problem.setVarType("head_"+startNode.getNodeID()+"_"+i, Double.class);
				problem.setVarLowerBound("head_"+startNode.getNodeID()+"_"+i, 0);
				
				problem.setVarType("k_"+startNode.getNodeID()+"_"+i, Boolean.class);
				
				problem.setVarType("yhead_"+startNode.getNodeID()+"_"+i, Double.class);
			}
			
			for(int j=0; j < pipeCost.size() ; j++){					
				problem.setVarType("l_"+i+"_"+j, Double.class);
				problem.setVarLowerBound("l_"+i+"_"+j, 0);
				
				problem.setVarType("y_"+i+"_"+j, Double.class);
				problem.setVarLowerBound("y_"+i+"_"+j, 0);
								
				if(pipe.isAllowParallel()){
					problem.setVarType("p_"+i+"_"+j, Boolean.class);
					problem.setVarType("yp_"+i+"_"+j, Boolean.class);
				}
			}
		}
				
		//s_i_j = 1 means ESR at ith node serves demand of jth node
		//d_i is the total demand served by ith esr
		//head_i is the head at ith node
		
		//e_i_j represents whether ith esr's cost is computed by jth esrcost table row
		//z_i_j = d_i * e_i_j
		
		// esr_i : height of esr at node i
		// ehead_i = elevation_i + esr_i - head_i
		for(int i : nodes.keySet()){			
			problem.setVarType("d_"+i, Double.class);
			
			problem.setVarType("esr_"+i, Double.class);
			problem.setVarLowerBound("esr_"+i, 0);
			
			problem.setVarType("besr_"+i, Boolean.class);
			
			problem.setVarType("head_"+i, Double.class);
			problem.setVarLowerBound("head_"+i, 0);
			
			problem.setVarType("ehead_"+i, Double.class);
			
			for(int j : nodes.keySet()){
				problem.setVarType("s_"+i+"_"+j, Boolean.class);
			}
						
			for(int j=0; j<esrCost.size();j++){
				problem.setVarType("e_"+i+"_"+j, Boolean.class);
				
				problem.setVarType("z_"+i+"_"+j, Double.class);
				problem.setVarLowerBound("z_"+i+"_"+j, 0);
			}
		}
	}
	
	//add variables to the ILP (includes zero demand nodes)
	private void addVariables_gen() throws Exception{
		// for each link i, for each commercial pipe j:
		// l_i_j : length of commercial pipe j of link i
		// p_i_j : boolean if ith link has jth pipe diameter in parallel
		// p_i : if ith link has no parallel pipe despite being allowed to have one
		//f_i denotes whether flow in pipe is primary = 1 or secondary = 0
		//y_i_j = f_i * l_i_j
		//yp_i_j = p_i_j * f_i
		
		// head_i is the incoming head at node i
		// introduce head_i_j : i is the start node, j is the id of the pipe
		// head_i_j is the head provided by i to the outgoing pipe j
		// if there is no ESR at i or if there is an ESR at i but water to pipe j still comes from the primary network, head_i_j = head_i
		// if there is an ESR at i and it is responsible for providing water to the pipe j, head_i_j = elevation_i + esr_i 
		// i.e. head provided is determined by the esr height and the location elevation
		// introduce head_i_j : is the start node, j is the id of the pipe
		// k_i_j = s_i_i && !f_j
		// head_i_j = k_i_j * (elevation_i + esr_i - head_i) + head_i
		// ehead_i = elevation_i + esr_i - head_i
		// yhead_i_j = k_i_j * ehead_i
		for(Pipe pipe : pipes.values()){
			int i = pipe.getPipeID();	
			problem.setVarType("f_"+i, Boolean.class);
			
			problem.setVarType("headloss_"+i, Double.class);
			problem.setVarLowerBound("headloss_"+i, 0);
			
			if(pipe.isAllowParallel()){
				problem.setVarType("p_"+i, Boolean.class);	
				problem.setVarType("yp_"+i, Boolean.class);
			}
			
			Node startNode = pipe.getStartNode();

			problem.setVarType("head_"+startNode.getNodeID()+"_"+i, Double.class);
			problem.setVarLowerBound("head_"+startNode.getNodeID()+"_"+i, 0);
			
			problem.setVarType("k_"+startNode.getNodeID()+"_"+i, Boolean.class);
			
			problem.setVarType("yhead_"+startNode.getNodeID()+"_"+i, Double.class);

			
			for(int j=0; j < pipeCost.size() ; j++){					
				problem.setVarType("l_"+i+"_"+j, Double.class);
				problem.setVarLowerBound("l_"+i+"_"+j, 0);
				
				problem.setVarType("y_"+i+"_"+j, Double.class);
				problem.setVarLowerBound("y_"+i+"_"+j, 0);
								
				if(pipe.isAllowParallel()){
					problem.setVarType("p_"+i+"_"+j, Boolean.class);
					problem.setVarType("yp_"+i+"_"+j, Boolean.class);
				}
			}
		}
				
		//s_i_j = 1 means ESR at ith node serves demand of jth node
		//d_i is the total demand served by ith esr
		//head_i is the head at ith node
		
		//e_i_j represents whether ith esr's cost is computed by jth esrcost table row
		//z_i_j = d_i * e_i_j
		
		// esr_i : height of esr at node i
		// ehead_i = elevation_i + esr_i - head_i
		for(int i : nodes.keySet()){			
			problem.setVarType("d_"+i, Double.class);
			
			problem.setVarType("esr_"+i, Double.class);
			problem.setVarLowerBound("esr_"+i, 0);
			
			problem.setVarType("besr_"+i, Boolean.class);
			
			problem.setVarType("head_"+i, Double.class);
			problem.setVarLowerBound("head_"+i, 0);
			
			problem.setVarType("ehead_"+i, Double.class);
			
			for(int j : nodes.keySet()){
				problem.setVarType("s_"+i+"_"+j, Boolean.class);
			}
						
			for(int j=0; j<esrCost.size();j++){
				problem.setVarType("e_"+i+"_"+j, Boolean.class);
				
				problem.setVarType("z_"+i+"_"+j, Double.class);
				problem.setVarLowerBound("z_"+i+"_"+j, 0);
			}
		}
	}
	
	//add variables to the ILP
	//considers l_i_j_0 and l_i_j_1 for secondary and primary networks)
	private void addVariables_gen3() throws Exception{
		// for each link i, for each commercial pipe j:
		// l_i_j_k : length of commercial pipe j of link i, k =1 if primary, 0 if secondary
		// p_i_j_k : boolean if ith link has jth pipe diameter in parallel, k =1 if primary, 0 if secondary
		// p_i : if ith link has no parallel pipe despite being allowed to have one
		//f_i denotes whether flow in pipe is primary = 1 or secondary = 0

		// head_i is the incoming head at node i
		// introduce head_i_j : i is the start node, j is the id of the pipe
		// head_i_j is the head provided by i to the outgoing pipe j
		// if there is no ESR at i or if there is an ESR at i but water to pipe j still comes from the primary network, head_i_j = head_i
		// if there is an ESR at i and it is responsible for providing water to the pipe j, head_i_j = elevation_i + esr_i 
		// i.e. head provided is determined by the esr height and the location elevation
		// introduce head_i_j : is the start node, j is the id of the pipe
		// k_i_j = s_i_i && !f_j
		// head_i_j = k_i_j * (elevation_i + esr_i - head_i) + head_i
		// ehead_i = elevation_i + esr_i - head_i
		// yhead_i_j = k_i_j * ehead_i
		for(Pipe pipe : pipes.values()){
			int i = pipe.getPipeID();	
			problem.setVarType("f_"+i, Boolean.class);
			
			problem.setVarType("headloss_"+i, Double.class);
			problem.setVarLowerBound("headloss_"+i, 0);
			
			if(pipe.isAllowParallel()){
				problem.setVarType("p_"+i+"_0", Boolean.class);
				problem.setVarType("p_"+i+"_1", Boolean.class);
			}
			
			Node startNode = pipe.getStartNode();

			problem.setVarType("head_"+startNode.getNodeID()+"_"+i, Double.class);
			problem.setVarLowerBound("head_"+startNode.getNodeID()+"_"+i, 0);
			
			problem.setVarType("k_"+startNode.getNodeID()+"_"+i, Boolean.class);
			
			problem.setVarType("yhead_"+startNode.getNodeID()+"_"+i, Double.class);

			
			for(int j=0; j < pipeCost.size() ; j++){					
				problem.setVarType("l_"+i+"_"+j+"_0", Double.class);
				problem.setVarLowerBound("l_"+i+"_"+j+"_0", 0);
				problem.setVarType("l_"+i+"_"+j+"_1", Double.class);
				problem.setVarLowerBound("l_"+i+"_"+j+"_1", 0);
								
				if(pipe.isAllowParallel()){
					problem.setVarType("p_"+i+"_"+j+"_0", Boolean.class);
					problem.setVarType("p_"+i+"_"+j+"_1", Boolean.class);
				}
			}
		}
				
		//s_i_j = 1 means ESR at ith node serves demand of jth node
		//d_i is the total demand served by ith esr
		//head_i is the head at ith node
		
		//e_i_j represents whether ith esr's cost is computed by jth esrcost table row
		//z_i_j = d_i * e_i_j
		
		// esr_i : height of esr at node i
		// ehead_i = elevation_i + esr_i - head_i
		for(int i : nodes.keySet()){			
			problem.setVarType("d_"+i, Double.class);
			
			problem.setVarType("esr_"+i, Double.class);
			problem.setVarLowerBound("esr_"+i, 0);
			
			problem.setVarType("besr_"+i, Boolean.class);
			
			problem.setVarType("head_"+i, Double.class);
			problem.setVarLowerBound("head_"+i, 0);
			
			problem.setVarType("ehead_"+i, Double.class);
			
			for(int j : nodes.keySet()){
				problem.setVarType("s_"+i+"_"+j, Boolean.class);
			}
						
			for(int j=0; j<esrCost.size();j++){
				problem.setVarType("e_"+i+"_"+j, Boolean.class);
				
				problem.setVarType("z_"+i+"_"+j, Double.class);
				problem.setVarLowerBound("z_"+i+"_"+j, 0);
			}
		}
	}
	
	//add variables to the ILP
	//includes pump related variables
	private void addVariables_gen4() throws Exception{
		// for each link i, for each commercial pipe j:
		// l_i_j_k : length of commercial pipe j of link i, k =1 if primary, 0 if secondary
		// p_i_j_k : boolean if ith link has jth pipe diameter in parallel, k =1 if primary, 0 if secondary
		// p_i : if ith link has no parallel pipe despite being allowed to have one
		// f_i denotes whether flow in pipe is primary = 1 or secondary = 0
		// pumphead_i : head provided by pump in pipe i
		// pumppower_i : power of pump installed in pipe i
		// ppumppower_i : primary power
		// spumppower_i : secondary power
		// f_pumphead_i : f_i * pumphead_i
		
		// head_i is the incoming head at node i
		// introduce head_i_j : i is the start node, j is the id of the pipe
		// head_i_j is the head provided by i to the outgoing pipe j
		// if there is no ESR at i or if there is an ESR at i but water to pipe j still comes from the primary network, head_i_j = head_i
		// if there is an ESR at i and it is responsible for providing water to the pipe j, head_i_j = elevation_i + esr_i 
		// i.e. head provided is determined by the esr height and the location elevation
		// introduce head_i_j : is the start node, j is the id of the pipe
		// k_i_j = s_i_i && !f_j
		// head_i_j = k_i_j * (elevation_i + esr_i - head_i) + head_i
		// ehead_i = elevation_i + esr_i - head_i
		// yhead_i_j = k_i_j * ehead_i
		for(Pipe pipe : pipes.values()){
			int i = pipe.getPipeID();	
			problem.setVarType("f_"+i, Boolean.class);
			
			problem.setVarType("headloss_"+i, Double.class);
			
			if(pumpGeneralProperties.pump_enabled){
				problem.setVarType("pumphead_"+i, Double.class);
				problem.setVarLowerBound("pumphead_"+i, 0);
				problem.setVarUpperBound("pumphead_"+i, maxPumpHead);
				
				
				problem.setVarType("pumppower_"+i, Double.class);
				problem.setVarLowerBound("pumppower_"+i, 0);
				
				problem.setVarType("ppumppower_"+i, Double.class);
				problem.setVarLowerBound("ppumppower_"+i, 0);
				
				problem.setVarType("spumppower_"+i, Double.class);
				problem.setVarLowerBound("spumppower_"+i, 0);
				
				problem.setVarType("f_pumphead_"+i, Double.class);
				problem.setVarLowerBound("f_pumphead_"+i, 0);
				
				problem.setVarType("pumphelper_"+i, Boolean.class);
			}
			
			if(pipe.isAllowParallel()){
				problem.setVarType("p_"+i+"_0", Boolean.class);
				problem.setVarType("p_"+i+"_1", Boolean.class);
			}
			
			Node startNode = pipe.getStartNode();

			problem.setVarType("head_"+startNode.getNodeID()+"_"+i, Double.class);
			problem.setVarLowerBound("head_"+startNode.getNodeID()+"_"+i, 0);
			
			problem.setVarType("k_"+startNode.getNodeID()+"_"+i, Boolean.class);
			
			problem.setVarType("yhead_"+startNode.getNodeID()+"_"+i, Double.class);

			
			for(int j=0; j < pipeCost.size() ; j++){					
				problem.setVarType("l_"+i+"_"+j+"_0", Double.class);
				problem.setVarLowerBound("l_"+i+"_"+j+"_0", 0);
				problem.setVarType("l_"+i+"_"+j+"_1", Double.class);
				problem.setVarLowerBound("l_"+i+"_"+j+"_1", 0);
								
				if(pipe.isAllowParallel()){
					problem.setVarType("p_"+i+"_"+j+"_0", Boolean.class);
					problem.setVarType("p_"+i+"_"+j+"_1", Boolean.class);
				}
			}
		}
		
		//s_i_j = 1 means ESR at ith node serves demand of jth node
		//d_i is the total demand served by ith esr
		//head_i is the head at ith node
		
		//e_i_j represents whether ith esr's cost is computed by jth esrcost table row
		//z_i_j = d_i * e_i_j
		
		// esr_i : height of esr at node i
		// ehead_i = elevation_i + esr_i - head_i
		for(int i : nodes.keySet()){			
			problem.setVarType("d_"+i, Double.class);
			
			problem.setVarType("esr_"+i, Double.class);
			problem.setVarLowerBound("esr_"+i, 0);
			
			problem.setVarType("besr_"+i, Boolean.class);
			
			problem.setVarType("head_"+i, Double.class);
			problem.setVarLowerBound("head_"+i, 0);
			
			problem.setVarType("ehead_"+i, Double.class);
			
			for(int j : nodes.keySet()){
				problem.setVarType("s_"+i+"_"+j, Boolean.class);
			}
						
			for(int j=0; j<esrCost.size();j++){
				problem.setVarType("e_"+i+"_"+j, Boolean.class);
				
				problem.setVarType("z_"+i+"_"+j, Double.class);
				problem.setVarLowerBound("z_"+i+"_"+j, 0);
			}
		}
	}
	
	//add variables to the ILP
	//does not include node variables s_i_j
	private void addVariables_gen5() throws Exception{
		// for each link i, for each commercial pipe j:
		// l_i_j_k : length of commercial pipe j of link i, k =1 if primary, 0 if secondary
		// p_i_j_k : boolean if ith link has jth pipe diameter in parallel, k =1 if primary, 0 if secondary
		// p_i : if ith link has no parallel pipe despite being allowed to have one
		// f_i denotes whether flow in pipe is primary = 1 or secondary = 0
		// pumphead_i : head provided by pump in pipe i
		// pumppower_i : power of pump installed in pipe i
		// ppumppower_i : primary power
		// spumppower_i : secondary power
		// f_pumphead_i : f_i * pumphead_i
		
		// head_i is the incoming head at node i
		// introduce head_i_j : i is the start node, j is the id of the pipe
		// head_i_j is the head provided by i to the outgoing pipe j
		// if there is no ESR at i or if there is an ESR at i but water to pipe j still comes from the primary network, head_i_j = head_i
		// if there is an ESR at i and it is responsible for providing water to the pipe j, head_i_j = elevation_i + esr_i 
		// i.e. head provided is determined by the esr height and the location elevation
		// introduce head_i_j : is the start node, j is the id of the pipe
		// k_i_j = s_i_i && !f_j
		// head_i_j = k_i_j * (elevation_i + esr_i - head_i) + head_i
		// ehead_i = elevation_i + esr_i - head_i
		// yhead_i_j = k_i_j * ehead_i
		for(Pipe pipe : pipes.values()){
			int i = pipe.getPipeID();	
			problem.setVarType("f_"+i, Boolean.class);
			
			problem.setVarType("headloss_"+i, Double.class);
			
			if(pumpGeneralProperties.pump_enabled){
				problem.setVarType("pumphead_"+i, Double.class);
				problem.setVarLowerBound("pumphead_"+i, 0);
				problem.setVarUpperBound("pumphead_"+i, maxPumpHead);
				
				
				problem.setVarType("pumppower_"+i, Double.class);
				problem.setVarLowerBound("pumppower_"+i, 0);
				
				problem.setVarType("ppumppower_"+i, Double.class);
				problem.setVarLowerBound("ppumppower_"+i, 0);
				
				problem.setVarType("spumppower_"+i, Double.class);
				problem.setVarLowerBound("spumppower_"+i, 0);
				
				problem.setVarType("f_pumphead_"+i, Double.class);
				problem.setVarLowerBound("f_pumphead_"+i, 0);
				
				problem.setVarType("pumphelper_"+i, Boolean.class);
			}
			
			if(pipe.isAllowParallel()){
				problem.setVarType("p_"+i+"_0", Boolean.class);
				problem.setVarType("p_"+i+"_1", Boolean.class);
			}
			
			Node startNode = pipe.getStartNode();

			problem.setVarType("head_"+startNode.getNodeID()+"_"+i, Double.class);
			problem.setVarLowerBound("head_"+startNode.getNodeID()+"_"+i, 0);
			
			problem.setVarType("k_"+startNode.getNodeID()+"_"+i, Boolean.class);
			
			problem.setVarType("yhead_"+startNode.getNodeID()+"_"+i, Double.class);

			
			for(int j=0; j < pipeCost.size() ; j++){					
				problem.setVarType("l_"+i+"_"+j+"_0", Double.class);
				problem.setVarLowerBound("l_"+i+"_"+j+"_0", 0);
				problem.setVarType("l_"+i+"_"+j+"_1", Double.class);
				problem.setVarLowerBound("l_"+i+"_"+j+"_1", 0);
								
				if(pipe.isAllowParallel()){
					problem.setVarType("p_"+i+"_"+j+"_0", Boolean.class);
					problem.setVarType("p_"+i+"_"+j+"_1", Boolean.class);
				}
			}
		}
		
		//d_i is the total demand served by ith esr
		//head_i is the head at ith node
		
		//e_i_j represents whether ith esr's cost is computed by jth esrcost table row
		//z_i_j = d_i * e_i_j
		
		// esr_i : height of esr at node i
		// ehead_i = elevation_i + esr_i - head_i
		for(int i : nodes.keySet()){			
			problem.setVarType("d_"+i, Double.class);
			
			problem.setVarType("esr_"+i, Double.class);
			problem.setVarLowerBound("esr_"+i, 0);
			
			problem.setVarType("besr_"+i, Boolean.class);
			
			problem.setVarType("head_"+i, Double.class);
			problem.setVarLowerBound("head_"+i, 0);
			
			problem.setVarType("ehead_"+i, Double.class);
					
			for(int j=0; j<esrCost.size();j++){
				problem.setVarType("e_"+i+"_"+j, Boolean.class);
				
				problem.setVarType("z_"+i+"_"+j, Double.class);
				problem.setVarLowerBound("z_"+i+"_"+j, 0);
			}
		}
	}
	
	//set the objective function of the ILP
	//includes ESR capital cost
	private void setObjectiveCost_esr() throws Exception{
		Linear linear = new Linear();
		int j=0;
		for(Pipe pipe : pipes.values()){
			int i = pipe.getPipeID();		
			j=0;
			for(PipeCost entry : pipeCost){				
				if(pipe.getDiameter()==0) // cost contributes only if diameter is to be computed
					linear.add(entry.getCost(), "l_"+i+"_"+j);
				if(pipe.isAllowParallel()){
					linear.add(entry.getCost() * pipe.getLength(), "p_"+i+"_"+j);
				}
				j++;
			}	
		}
		
		for(Node node : nodes.values()){
			int i = node.getNodeID();
			if(node.getDemand()>0){
				j=0;
				for(EsrCost esr : esrCost){
					linear.add(esr.getBaseCost() - esr.getMinCapacity()*esr.getUnitCost(),"e_"+i+"_"+j);
					linear.add(esr.getUnitCost(),"z_"+i+"_"+j);
					j++;
				}
			}
		}
		
		//problem.setObjective(linear, OptType.MIN);
		problem.setObjective(linear, false);
	}
	
	//set the objective function of the ILP
	//allows nodes with zero demands to have ESRs 
	private void setObjectiveCost_esr_gen() throws Exception{
		Linear linear = new Linear();
		int j=0;
		for(Pipe pipe : pipes.values()){
			int i = pipe.getPipeID();		
			j=0;
			for(PipeCost entry : pipeCost){				
				if(pipe.getDiameter()==0) // cost contributes only if diameter is to be computed
					linear.add(entry.getCost(), "l_"+i+"_"+j);
				if(pipe.isAllowParallel()){
					linear.add(entry.getCost() * pipe.getLength(), "p_"+i+"_"+j);
				}
				j++;
			}	
		}
		
		for(Node node : nodes.values()){
			int i = node.getNodeID();
				j=0;
				for(EsrCost esr : esrCost){
					linear.add(esr.getBaseCost() - esr.getMinCapacity()*esr.getUnitCost(),"e_"+i+"_"+j);
					linear.add(esr.getUnitCost(),"z_"+i+"_"+j);
					j++;
				}
		}
		
		//esr height factor
		double minesrcost = getCost(totalDemand, esrCost);
		double nodenumber = nodes.size();
		double avgesrcost = 0.01*minesrcost/nodenumber;
		for(Node node : nodes.values()){
			int i = node.getNodeID();
			linear.add(avgesrcost,"esr_"+i);
		}
		
		//problem.setObjective(linear, OptType.MIN);
		problem.setObjective(linear, false);
	}
	
	//set the objective function of the ILP
	//updated to consider l_i_j_k variables
	private void setObjectiveCost_esr_gen3() throws Exception{
		Linear linear = new Linear();
		int j=0;
		for(Pipe pipe : pipes.values()){
			int i = pipe.getPipeID();		
			j=0;
			for(PipeCost entry : pipeCost){				
				if(pipe.getDiameter()==0){ // cost contributes only if diameter is to be computed
					linear.add(entry.getCost(), "l_"+i+"_"+j+"_0");
					linear.add(entry.getCost(), "l_"+i+"_"+j+"_1");
				}
				if(pipe.isAllowParallel()){
					linear.add(entry.getCost() * pipe.getLength(), "p_"+i+"_"+j+"_0");
					linear.add(entry.getCost() * pipe.getLength(), "p_"+i+"_"+j+"_1");
				}
				j++;
			}	
		}
		
		for(Node node : nodes.values()){
			int i = node.getNodeID();
				j=0;
				for(EsrCost esr : esrCost){
					linear.add(esr.getBaseCost() - esr.getMinCapacity()*esr.getUnitCost(),"e_"+i+"_"+j);
					linear.add(esr.getUnitCost(),"z_"+i+"_"+j);
					j++;
				}
		}
		
		
		//esr height factor
		double minesrcost = getCost(totalDemand, esrCost);
		double nodenumber = nodes.size();
		double avgesrcost = 0.01*minesrcost/nodenumber;
		for(Node node : nodes.values()){
			int i = node.getNodeID();
			linear.add(avgesrcost,"esr_"+i);
		}
		
		//problem.setObjective(linear, OptType.MIN);
		problem.setObjective(linear, false);
	}
	
	//set the objective function of the ILP
	//includes pump cost
	private void setObjectiveCost_esr_gen4() throws Exception{
		Linear linear = new Linear();
		int j=0;
		for(Pipe pipe : pipes.values()){
			int i = pipe.getPipeID();		
			j=0;
			for(PipeCost entry : pipeCost){				
				if(pipe.getDiameter()==0){ // cost contributes only if diameter is to be computed
					linear.add(entry.getCost(), "l_"+i+"_"+j+"_0");
					linear.add(entry.getCost(), "l_"+i+"_"+j+"_1");
				}
				if(pipe.isAllowParallel()){
					linear.add(entry.getCost() * pipe.getLength(), "p_"+i+"_"+j+"_0");
					linear.add(entry.getCost() * pipe.getLength(), "p_"+i+"_"+j+"_1");
				}
				j++;
			}
			
			// energy cost = power * hours
			// power contains flow term
			// flow = base_flow *base_hours/hours
			// hours term cancel
			// therefore for both primary and secondary pipes
			// energy cost = power as per base flow * base_hours
			
			// flow_m3h = flow*3.6;
			//power = density of water (1000) * g (9.80665) * flow (in m3h = flow*3.6) / (3.6*10^6) * efficiency
			// power = 9.81 * flow (in lps) / (1000 * efficiency)
			// energy cost = power * no of hours * cost per kw
			if(pumpGeneralProperties.pump_enabled){
				double presentvaluefactor = Util.presentValueFactor(pumpGeneralProperties.discount_rate, pumpGeneralProperties.inflation_rate, pumpGeneralProperties.design_lifetime);
				//double pumpcoeffecient = 9.80665*365*presentvaluefactor*pipe.getFlow()*generalProperties.supply_hours*pumpGeneralProperties.energycost_per_kw*pumpGeneralProperties.energycost_factor/(1000*pumpGeneralProperties.efficiency);
				double primarycoefficient = presentvaluefactor*365*generalProperties.supply_hours*pumpGeneralProperties.energycost_per_kwh;
				double secondarycoefficient = presentvaluefactor*365*esrGeneralProperties.secondary_supply_hours*pumpGeneralProperties.energycost_per_kwh;
				
				// energy cost
				linear.add(primarycoefficient, "ppumppower_"+i);
				linear.add(secondarycoefficient, "spumppower_"+i);
				
				//capital cost
				linear.add(pumpGeneralProperties.capitalcost_per_kw, "pumppower_"+i);
			}
		}
		
		for(Node node : nodes.values()){
			int i = node.getNodeID();
				j=0;
				for(EsrCost esr : esrCost){
					linear.add(esr.getBaseCost() - esr.getMinCapacity()*esr.getUnitCost(),"e_"+i+"_"+j);
					linear.add(esr.getUnitCost(),"z_"+i+"_"+j);
					j++;
				}
		}
		
		
		//esr height factor
		double minesrcost = getCost(totalDemand, esrCost);
		double nodenumber = nodes.size();
		double avgesrcost = 0.01*minesrcost/nodenumber;
		for(Node node : nodes.values()){
			int i = node.getNodeID();
			linear.add(avgesrcost,"esr_"+i);
		}
		
		//problem.setObjective(linear, OptType.MIN);
		problem.setObjective(linear, false);
	}
	
	//set the objective function of the ILP
	//removal of s_i_j variables
	private void setObjectiveCost_esr_gen5() throws Exception{
		Linear linear = new Linear();
		int j=0;
		for(Pipe pipe : pipes.values()){
			int i = pipe.getPipeID();		
			j=0;
			for(PipeCost entry : pipeCost){				
				if(pipe.getDiameter()==0){ // cost contributes only if diameter is to be computed
					linear.add(entry.getCost(), "l_"+i+"_"+j+"_0");
					linear.add(entry.getCost(), "l_"+i+"_"+j+"_1");
				}
				if(pipe.isAllowParallel()){
					linear.add(entry.getCost() * pipe.getLength(), "p_"+i+"_"+j+"_0");
					linear.add(entry.getCost() * pipe.getLength(), "p_"+i+"_"+j+"_1");
				}
				j++;
			}
			
			// energy cost = power * hours
			// power contains flow term
			// flow = base_flow *base_hours/hours
			// hours term cancel
			// therefore for both primary and secondary pipes
			// energy cost = power as per base flow * base_hours
			
			// flow_m3h = flow*3.6;
			//power = density of water (1000) * g (9.80665) * flow (in m3h = flow*3.6) / (3.6*10^6) * efficiency
			// power = 9.81 * flow (in lps) / (1000 * efficiency)
			// energy cost = power * no of hours * cost per kw
			if(pumpGeneralProperties.pump_enabled){
				double presentvaluefactor = Util.presentValueFactor(pumpGeneralProperties.discount_rate, pumpGeneralProperties.inflation_rate, pumpGeneralProperties.design_lifetime);
				//double pumpcoeffecient = 9.80665*365*presentvaluefactor*pipe.getFlow()*generalProperties.supply_hours*pumpGeneralProperties.energycost_per_kw*pumpGeneralProperties.energycost_factor/(1000*pumpGeneralProperties.efficiency);
				double primarycoefficient = presentvaluefactor*365*generalProperties.supply_hours*pumpGeneralProperties.energycost_per_kwh;
				double secondarycoefficient = presentvaluefactor*365*esrGeneralProperties.secondary_supply_hours*pumpGeneralProperties.energycost_per_kwh;
				
				// energy cost
				linear.add(primarycoefficient, "ppumppower_"+i);
				linear.add(secondarycoefficient, "spumppower_"+i);
				
				//capital cost
				linear.add(pumpGeneralProperties.capitalcost_per_kw, "pumppower_"+i);
			}
		}
		
		for(Node node : nodes.values()){
			int i = node.getNodeID();
				j=0;
				for(EsrCost esr : esrCost){
					linear.add(esr.getBaseCost() - esr.getMinCapacity()*esr.getUnitCost(),"e_"+i+"_"+j);
					linear.add(esr.getUnitCost(),"z_"+i+"_"+j);
					j++;
				}
		}
		
		
		//esr height factor
		double minesrcost = getCost(totalDemand, esrCost);
		double nodenumber = nodes.size();
		double avgesrcost = 0.01*minesrcost/nodenumber;
		for(Node node : nodes.values()){
			int i = node.getNodeID();
			linear.add(avgesrcost,"esr_"+i);
		}
		
		//problem.setObjective(linear, OptType.MIN);
		problem.setObjective(linear, false);
	}
	
	//add constraints for the pipes in the network
	//includes constraints related to headloss/ water speed limits and assignment of commercial pipe diameters
	private void setPipeConstraints() throws Exception{
		for(Pipe pipe : pipes.values()){
			int i = pipe.getPipeID();
			Linear linear;
			boolean atleastOneHeadlossCorrect = false;
			boolean atleastOneSpeedCorrect = false;
			
			//true when pipe has already been assigned a diameter by the user
			if(pipe.getDiameter()!=0){
				int j=0;
				boolean diameterFromPipeCost = false;
				for(PipeCost entry : pipeCost){
					linear = new Linear();
					linear.add(1, "l_"+i+"_"+j);
					if(pipe.getDiameter() == entry.getDiameter()){
						problem.add(new Constraint("", linear, "=", pipe.getLength()));
						diameterFromPipeCost = true;
					}
					else
						problem.add(new Constraint("", linear, "=", 0));				
					j++;
				}
				
				//if diameter provided by user is not defined in the pipecost table
				if(!diameterFromPipeCost)
					throw new Exception("The custom diameter: " + pipe.getDiameter() + " with roughness: "+ pipe.getRoughness() +" does not belong to commercial pipes"); 
				
				//if pipe is allowed to have a parallel pipe, add corresponding constraints
				if(pipe.isAllowParallel()){
					linear = new Linear();			
					for(int j1=0;j1<pipeCost.size();j1++){
						linear.add(1, "p_"+i+"_"+j1);
					}
					linear.add(1, "p_"+i);
					problem.add(new Constraint("", linear, "=", 1));	
					
					j=0;
					for(PipeCost entry : pipeCost){
						double flow = pipe.getFlow();
									
						//when parallel pipe of dia j, what is flow in primary
						flow = flow / (1 + (entry.getRoughness()/pipe.getRoughness())*Math.pow(entry.getDiameter()/pipe.getDiameter(), 4.87/1.852));						
						double headloss = Util.HWheadLoss(flow, pipe.getRoughness(), pipe.getDiameter());
						
						if(headloss < generalProperties.min_hl_perkm/1000 || headloss > generalProperties.max_hl_perkm/1000){
							linear = new Linear();
							linear.add(1, "p_"+i+"_"+j);
							problem.add(new Constraint("", linear, "=", 0));
							//System.out.println(i+" "+entry.getDiameter()+" "+headloss);
						}
						else
							atleastOneHeadlossCorrect = true;
						
						double parallel_flow = pipe.getFlow() - flow;
						double main_speed = Util.waterSpeed(flow, pipe.getDiameter()); 
						double parallel_speed = Util.waterSpeed(parallel_flow, entry.getDiameter()); 
						if( generalProperties.max_water_speed > 0 && (parallel_speed > generalProperties.max_water_speed || main_speed > generalProperties.max_water_speed)){
							linear = new Linear();
							linear.add(1, "p_"+i+"_"+j);
							problem.add(new Constraint("", linear, "=", 0));
						}
						else
							atleastOneSpeedCorrect = true;
						j++;
					}
					
					double flow = pipe.getFlow();
					double headloss = Util.HWheadLoss(flow, pipe.getRoughness(), pipe.getDiameter());
					
					if(headloss < generalProperties.min_hl_perkm/1000 || headloss > generalProperties.max_hl_perkm/1000){
						linear = new Linear();
						linear.add(1, "p_"+i);
						problem.add(new Constraint("", linear, "=", 0));
						//System.out.println(i+" "+pipe.getDiameter()+" "+headloss);
					}
					else
						atleastOneHeadlossCorrect = true;
					
					double speed = Util.waterSpeed(flow, pipe.getDiameter());
					if(generalProperties.max_water_speed > 0 && speed > generalProperties.max_water_speed){
						linear = new Linear();
						linear.add(1, "p_"+i);
						problem.add(new Constraint("", linear, "=", 0));
					}
					else
						atleastOneSpeedCorrect = true;
				}
				else{
					double flow = pipe.getFlow();
					double headloss = Util.HWheadLoss(flow, pipe.getRoughness(), pipe.getDiameter());
					
					if(headloss < generalProperties.min_hl_perkm/1000 || headloss > generalProperties.max_hl_perkm/1000){
						
					}
					else
						atleastOneHeadlossCorrect = true;
										
					double primary_speed = Util.waterSpeed(flow, pipe.getDiameter());
					
					if(generalProperties.max_water_speed > 0 && primary_speed > generalProperties.max_water_speed){
						
					}
					else
						atleastOneSpeedCorrect = true;					
				}
			}
			else{
				linear = new Linear();
				for(int j=0;j<pipeCost.size();j++){
					linear.add(1, "l_"+i+"_"+j);
				}
				problem.add(new Constraint("", linear, "=", pipe.getLength()));	
				
				int j=0;
				
				for(PipeCost entry : pipeCost){
					double flow = pipe.getFlow();
					double headloss = Util.HWheadLoss(flow, entry.getRoughness(), entry.getDiameter());
					double speed = Util.waterSpeed(flow, entry.getDiameter()); 
					
					if(headloss < generalProperties.min_hl_perkm/1000 || headloss > generalProperties.max_hl_perkm/1000){
						linear = new Linear();
						linear.add(1, "l_"+i+"_"+j);
						problem.add(new Constraint("", linear, "=", 0));
						//System.out.println(i+" "+entry.getDiameter()+" "+headloss);
					}
					else
						atleastOneHeadlossCorrect = true;
					
					if(generalProperties.max_water_speed > 0 && speed > generalProperties.max_water_speed){
						linear = new Linear();
						linear.add(1, "l_"+i+"_"+j);
						problem.add(new Constraint("", linear, "=", 0));
					}
					else
						atleastOneSpeedCorrect = true;
					j++;
				}
			}
			
			//if none of the commercial pipe diameters can satisfy headloss/speed constraints for the pipe, terminate optimization and inform user
			if(!atleastOneHeadlossCorrect)
				throw new Exception("Headloss for pipe: " + pipe.getPipeID() + " cannot be set within the min/max headloss values.");
			
			if(!atleastOneSpeedCorrect)
				throw new Exception("Speed of water in pipe: " + pipe.getPipeID() + " cannot be set below the max value.");
		}
		//System.out.println(problem);
	}
	
	//add constraints for the pipes in the network
	//includes constraints related to headloss/ water speed limits and assignment of commercial pipe diameters
	//includes consideration for secondary network, if ESR optimization is enabled	
	private void setPipeConstraints_model3() throws Exception{
		for(Pipe pipe : pipes.values()){
			int i = pipe.getPipeID();
			Linear linear;
			boolean atleastOneHeadlossCorrect = false;
			boolean atleastOneSpeedCorrect = false;
			
			//true when pipe has already been assigned a diameter by the user
			if(pipe.getDiameter()!=0){
				int j=0;
				int fixedj = 0;
				boolean diameterFromPipeCost = false;
				for(PipeCost entry : pipeCost){
					linear = new Linear();
					linear.add(1, "l_"+i+"_"+j);
					if(pipe.getDiameter() == entry.getDiameter()){
						problem.add(new Constraint("", linear, "=", pipe.getLength()));
						diameterFromPipeCost = true;
						fixedj = j;
					}
					else
						problem.add(new Constraint("", linear, "=", 0));				
					j++;
				}
				
				//if diameter provided by user is not defined in the pipecost table
				if(!diameterFromPipeCost)
					throw new Exception("The custom diameter: " + pipe.getDiameter() + " with roughness: "+ pipe.getRoughness() +" does not belong to commercial pipes"); 
				
				//if pipe is allowed to have a parallel pipe, add corresponding constraints
				if(pipe.isAllowParallel()){
					linear = new Linear();			
					for(int j1=0;j1<pipeCost.size();j1++){
						linear.add(1, "p_"+i+"_"+j1);
					}
					linear.add(1, "p_"+i);
					problem.add(new Constraint("", linear, "=", 1));	
					
					j=0;
					for(PipeCost entry : pipeCost){
						double flow = pipe.getFlow();
						
						//when parallel pipe of dia j, what is flow in primary
						flow = flow / (1 + (entry.getRoughness()/pipe.getRoughness())*Math.pow(entry.getDiameter()/pipe.getDiameter(), 4.87/1.852));						
						double headloss = Util.HWheadLoss(flow, pipe.getRoughness(), pipe.getDiameter());
						double secondaryHeadloss = Util.HWheadLoss(flow * secondaryFlowFactor, pipe.getRoughness(), pipe.getDiameter());
												
						if(headloss < generalProperties.min_hl_perkm/1000 || headloss > generalProperties.max_hl_perkm/1000){
							linear = new Linear();
							linear.add(1, "yp_"+i+"_"+j);
							problem.add(new Constraint("", linear, "=", 0));
							//System.out.println(i+" "+entry.getDiameter()+" "+headloss);
						}
						else
							atleastOneHeadlossCorrect = true;
						
						if(secondaryHeadloss < generalProperties.min_hl_perkm/1000 || secondaryHeadloss > generalProperties.max_hl_perkm/1000){
							linear = new Linear();
							linear.add(1, "yp_"+i+"_"+j);
							linear.add(-1, "p_"+i+"_"+j);
							problem.add(new Constraint("", linear, "=", 0));
							//System.out.println(i+" "+entry.getDiameter()+" "+headloss);
						}
						else
							atleastOneHeadlossCorrect = true;
						
						double parallel_flow = pipe.getFlow() - flow;
						double main_speed = Util.waterSpeed(flow, pipe.getDiameter()); 
						double parallel_speed = Util.waterSpeed(parallel_flow, entry.getDiameter()); 
						
						double main_secondary_speed = Util.waterSpeed(flow*secondaryFlowFactor, pipe.getDiameter()); 
						double parallel_secondary_speed = Util.waterSpeed(parallel_flow*secondaryFlowFactor, entry.getDiameter());
						
						
						if(generalProperties.max_water_speed > 0 && (parallel_speed > generalProperties.max_water_speed || main_speed > generalProperties.max_water_speed)){
							linear = new Linear();
							linear.add(1, "yp_"+i+"_"+j);
							problem.add(new Constraint("", linear, "=", 0));
						}
						else
							atleastOneSpeedCorrect = true;
						
						if(generalProperties.max_water_speed > 0 && (parallel_secondary_speed > generalProperties.max_water_speed || main_secondary_speed > generalProperties.max_water_speed)){
							linear = new Linear();
							linear.add(1, "yp_"+i+"_"+j);
							linear.add(-1, "p_"+i+"_"+j);
							problem.add(new Constraint("", linear, "=", 0));
						}
						else
							atleastOneSpeedCorrect = true;
						
						j++;
					}
					
					double flow = pipe.getFlow();
					double headloss = Util.HWheadLoss(flow, pipe.getRoughness(), pipe.getDiameter());
					double secondaryHeadloss = Util.HWheadLoss(flow * secondaryFlowFactor, pipe.getRoughness(), pipe.getDiameter());
					
					if(headloss < generalProperties.min_hl_perkm/1000 || headloss > generalProperties.max_hl_perkm/1000){
						linear = new Linear();
						linear.add(1, "yp_"+i);
						problem.add(new Constraint("", linear, "=", 0));
						//System.out.println(i+" "+pipe.getDiameter()+" "+headloss);
					}
					else
						atleastOneHeadlossCorrect = true;
					
					if(secondaryHeadloss < generalProperties.min_hl_perkm/1000 || secondaryHeadloss > generalProperties.max_hl_perkm/1000){
						linear = new Linear();
						linear.add(1, "yp_"+i);
						linear.add(-1, "p_"+i);
						problem.add(new Constraint("", linear, "=", 0));
						//System.out.println(i+" "+pipe.getDiameter()+" "+headloss);
					}
					else
						atleastOneHeadlossCorrect = true;
					
					double speed = Util.waterSpeed(flow, pipe.getDiameter());
					double secondary_speed = Util.waterSpeed(flow*secondaryFlowFactor, pipe.getDiameter());
					
					if(generalProperties.max_water_speed > 0 && speed > generalProperties.max_water_speed){
						linear = new Linear();
						linear.add(1, "yp_"+i);
						problem.add(new Constraint("", linear, "=", 0));
					}
					else
						atleastOneSpeedCorrect = true;
					
					if(generalProperties.max_water_speed > 0 && secondary_speed > generalProperties.max_water_speed){
						linear = new Linear();
						linear.add(1, "yp_"+i);
						linear.add(-1, "p_"+i);
						problem.add(new Constraint("", linear, "=", 0));
					}
					else
						atleastOneSpeedCorrect = true;
				}
				else{
					double flow = pipe.getFlow();
					double headloss = Util.HWheadLoss(flow, pipe.getRoughness(), pipe.getDiameter());
					double secondaryHeadloss = Util.HWheadLoss(flow * secondaryFlowFactor, pipe.getRoughness(), pipe.getDiameter());
					
					if(headloss < generalProperties.min_hl_perkm/1000 || headloss > generalProperties.max_hl_perkm/1000){
						linear = new Linear();
						linear.add(1, "y_"+i+"_"+fixedj);
						problem.add(new Constraint("", linear, "=", 0));
					}
					else
						atleastOneHeadlossCorrect = true;
					
					if(secondaryHeadloss < generalProperties.min_hl_perkm/1000 || secondaryHeadloss > generalProperties.max_hl_perkm/1000){
						linear = new Linear();
						linear.add(1, "y_"+i+"_"+fixedj);
						linear.add(-1, "l_"+i+"_"+fixedj);
						problem.add(new Constraint("", linear, "=", 0));
					}
					else
						atleastOneHeadlossCorrect = true;
										
					double primary_speed = Util.waterSpeed(flow, pipe.getDiameter());
					double secondary_speed = Util.waterSpeed(flow * secondaryFlowFactor, pipe.getDiameter());
					
					
					if(generalProperties.max_water_speed > 0 && primary_speed > generalProperties.max_water_speed){
						linear = new Linear();
						linear.add(1, "y_"+i+"_"+fixedj);
						problem.add(new Constraint("", linear, "=", 0));
					}
					else
						atleastOneSpeedCorrect = true;	
					
					if(generalProperties.max_water_speed > 0 && secondary_speed > generalProperties.max_water_speed){
						linear = new Linear();
						linear.add(1, "y_"+i+"_"+fixedj);
						linear.add(-1, "l_"+i+"_"+fixedj);
						problem.add(new Constraint("", linear, "=", 0));
					}
					else
						atleastOneSpeedCorrect = true;
				}
			}
			else{
				linear = new Linear();
				for(int j=0;j<pipeCost.size();j++){
					linear.add(1, "l_"+i+"_"+j);
				}
				problem.add(new Constraint("", linear, "=", pipe.getLength()));	
				
				int j=0;
				
				for(PipeCost entry : pipeCost){
					double flow = pipe.getFlow();
					double headloss = Util.HWheadLoss(flow, entry.getRoughness(), entry.getDiameter());
					double secondaryHeadloss = Util.HWheadLoss(flow * secondaryFlowFactor, entry.getRoughness(), entry.getDiameter());
					double speed = Util.waterSpeed(flow, entry.getDiameter()); 
					double secondary_speed = Util.waterSpeed(flow * secondaryFlowFactor, entry.getDiameter()); 
							
					if(headloss < generalProperties.min_hl_perkm/1000 || headloss > generalProperties.max_hl_perkm/1000){
						linear = new Linear();
						linear.add(1, "y_"+i+"_"+j);
						problem.add(new Constraint("", linear, "=", 0));
						//System.out.println(i+" "+entry.getDiameter()+" "+headloss);
					}
					else
						atleastOneHeadlossCorrect = true;
					
					if(secondaryHeadloss < generalProperties.min_hl_perkm/1000 || secondaryHeadloss > generalProperties.max_hl_perkm/1000){
						linear = new Linear();
						linear.add(1, "y_"+i+"_"+j);
						linear.add(-1, "l_"+i+"_"+j);
						problem.add(new Constraint("", linear, "=", 0));
						//System.out.println(i+" "+entry.getDiameter()+" "+headloss);
					}
					else
						atleastOneHeadlossCorrect = true;
					
					if(generalProperties.max_water_speed > 0 && speed > generalProperties.max_water_speed){
						linear = new Linear();
						linear.add(1, "y_"+i+"_"+j);
						problem.add(new Constraint("", linear, "=", 0));
					}
					else
						atleastOneSpeedCorrect = true;
					
					if(generalProperties.max_water_speed > 0 && secondary_speed > generalProperties.max_water_speed){
						linear = new Linear();
						linear.add(1, "y_"+i+"_"+j);
						linear.add(-1, "l_"+i+"_"+j);
						problem.add(new Constraint("", linear, "=", 0));
					}
					else
						atleastOneSpeedCorrect = true;
					
					j++;
				}
			}
			
			//if none of the commercial pipe diameters can satisfy headloss/speed constraints for the pipe, terminate optimization and inform user
			if(!atleastOneHeadlossCorrect)
				throw new Exception("Headloss for pipe: " + pipe.getPipeID() + " cannot be set within the min/max headloss values.");
			
			if(!atleastOneSpeedCorrect)
				throw new Exception("Speed of water in pipe: " + pipe.getPipeID() + " cannot be set below the max value.");
		}
		//System.out.println(problem);
	}
	
	//add constraints for the pipes in the network
	//includes constraints related to headloss/ water speed limits and assignment of commercial pipe diameters
	//updated with l_i_j_k instead of l_i_j
	private void setPipeConstraints_gen3() throws Exception{		
		for(Pipe pipe : pipes.values()){
			int i = pipe.getPipeID();
			Linear linear;
			boolean atleastOneHeadlossCorrect = false;
			boolean atleastOneSpeedCorrect = false;
			
			linear = new Linear();
			Linear linear2 = new Linear();
			for(int j=0;j<pipeCost.size();j++){
				linear.add(1, "l_"+i+"_"+j+"_0");
				linear2.add(1, "l_"+i+"_"+j+"_1");
			}
			linear.add(pipe.getLength(), "f_"+i);
			linear2.add(-1*pipe.getLength(), "f_"+i);
			problem.add(new Constraint("", linear, "=", pipe.getLength()));	
			problem.add(new Constraint("", linear2, "=", 0));
			
			if(pipe.getDiameter()!=0){
				int j=0;
				int fixedj = 0;
				boolean diameterFromPipeCost = false;
				for(PipeCost entry : pipeCost){
					linear = new Linear();
					linear.add(1, "l_"+i+"_"+j+"_0");
					linear.add(1, "l_"+i+"_"+j+"_1");
					if(pipe.getDiameter() == entry.getDiameter()){
						problem.add(new Constraint("", linear, "=", pipe.getLength()));
						diameterFromPipeCost = true;
						fixedj = j;
					}
					else
						problem.add(new Constraint("", linear, "=", 0));				
					j++;
				}
				if(!diameterFromPipeCost)
					throw new Exception("The custom diameter: " + pipe.getDiameter() + " with roughness: "+ pipe.getRoughness() +" does not belong to commercial pipes"); 
				
				if(pipe.isAllowParallel()){
//					linear = new Linear();			
//					for(int j1=0;j1<pipeCost.size();j1++){
//						linear.add(1, "p_"+i+"_"+j1+"_0");
//						linear.add(1, "p_"+i+"_"+j1+"_1");
//					}
//					linear.add(1, "p_"+i+"_0");
//					linear.add(1, "p_"+i+"_1");
//					problem.add(new Constraint("", linear, "=", 1));	
					
					linear = new Linear();
					linear2 = new Linear();
					for(int j1=0;j1<pipeCost.size();j1++){
						linear.add(1, "p_"+i+"_"+j1+"_0");
						linear2.add(1, "p_"+i+"_"+j1+"_1");
					}
					linear.add(1, "p_"+i+"_0");
					linear2.add(1, "p_"+i+"_1");
					linear.add(1, "f_"+i);
					linear2.add(-1, "f_"+i);
					problem.add(new Constraint("", linear, "=", 1));
					problem.add(new Constraint("", linear2, "=", 0));
					
					j=0;
					for(PipeCost entry : pipeCost){
						double flow = pipe.getFlow();
									
						//when parallel pipe of dia j, what is flow in primary
						flow = flow / (1 + (entry.getRoughness()/pipe.getRoughness())*Math.pow(entry.getDiameter()/pipe.getDiameter(), 4.87/1.852));						
						double headloss = Util.HWheadLoss(flow, pipe.getRoughness(), pipe.getDiameter());
						double secondaryHeadloss = Util.HWheadLoss(flow * secondaryFlowFactor, pipe.getRoughness(), pipe.getDiameter());
						
						if(headloss < generalProperties.min_hl_perkm/1000 || headloss > generalProperties.max_hl_perkm/1000){
							linear = new Linear();
							linear.add(1, "p_"+i+"_"+j+"_1");
							problem.add(new Constraint("", linear, "=", 0));
						}
						else
							atleastOneHeadlossCorrect = true;
						
						if(secondaryHeadloss < generalProperties.min_hl_perkm/1000 || secondaryHeadloss > generalProperties.max_hl_perkm/1000){
							linear = new Linear();
							linear.add(1, "p_"+i+"_"+j+"_0");
							problem.add(new Constraint("", linear, "=", 0));
						}
						else
							atleastOneHeadlossCorrect = true;
						
						double parallel_flow = pipe.getFlow() - flow;
						double main_primary_speed = Util.waterSpeed(flow, pipe.getDiameter()); 
						double parallel_primary_speed = Util.waterSpeed(parallel_flow, entry.getDiameter());
						
						double main_secondary_speed = Util.waterSpeed(flow*secondaryFlowFactor, pipe.getDiameter()); 
						double parallel_secondary_speed = Util.waterSpeed(parallel_flow*secondaryFlowFactor, entry.getDiameter());
						
						
						if( generalProperties.max_water_speed > 0 && (main_primary_speed > generalProperties.max_water_speed || parallel_primary_speed > generalProperties.max_water_speed)){
							linear = new Linear();
							linear.add(1, "p_"+i+"_"+j+"_1");
							problem.add(new Constraint("", linear, "=", 0));
						}
						else
							atleastOneSpeedCorrect = true;
						
						if( generalProperties.max_water_speed > 0 && (main_secondary_speed > generalProperties.max_water_speed || parallel_secondary_speed > generalProperties.max_water_speed)){
							linear = new Linear();
							linear.add(1, "p_"+i+"_"+j+"_0");
							problem.add(new Constraint("", linear, "=", 0));
						}
						else
							atleastOneSpeedCorrect = true;
						
						j++;
					}
					
					double flow = pipe.getFlow();
					double headloss = Util.HWheadLoss(flow, pipe.getRoughness(), pipe.getDiameter());
					double secondaryHeadloss = Util.HWheadLoss(flow * secondaryFlowFactor, pipe.getRoughness(), pipe.getDiameter());
					
					if(headloss < generalProperties.min_hl_perkm/1000 || headloss > generalProperties.max_hl_perkm/1000){
						linear = new Linear();
						linear.add(1, "p_"+i+"_1");
						problem.add(new Constraint("", linear, "=", 0));
					}
					else
						atleastOneHeadlossCorrect = true;
					
					if(secondaryHeadloss < generalProperties.min_hl_perkm/1000 || secondaryHeadloss > generalProperties.max_hl_perkm/1000){
						linear = new Linear();
						linear.add(1, "p_"+i+"_0");
						problem.add(new Constraint("", linear, "=", 0));
					}
					else
						atleastOneHeadlossCorrect = true;
					
					double primary_speed = Util.waterSpeed(flow, pipe.getDiameter());
					double secondary_speed = Util.waterSpeed(flow*secondaryFlowFactor, pipe.getDiameter());
					
					if(generalProperties.max_water_speed > 0 && primary_speed > generalProperties.max_water_speed){
						linear = new Linear();
						linear.add(1, "p_"+i+"_1");
						problem.add(new Constraint("", linear, "=", 0));
					}
					else
						atleastOneSpeedCorrect = true;
					
					if(generalProperties.max_water_speed > 0 && secondary_speed > generalProperties.max_water_speed){
						linear = new Linear();
						linear.add(1, "p_"+i+"_0");
						problem.add(new Constraint("", linear, "=", 0));
					}
					else
						atleastOneSpeedCorrect = true;
					
				}
				else{
					double flow = pipe.getFlow();
					double headloss = Util.HWheadLoss(flow, pipe.getRoughness(), pipe.getDiameter());
					double secondaryHeadloss = Util.HWheadLoss(flow * secondaryFlowFactor, pipe.getRoughness(), pipe.getDiameter());
					
					if(headloss < generalProperties.min_hl_perkm/1000 || headloss > generalProperties.max_hl_perkm/1000){
						linear = new Linear();
						linear.add(1, "l_"+i+"_"+fixedj+"_1");
						problem.add(new Constraint("", linear, "=", 0));
					}
					else
						atleastOneHeadlossCorrect = true;
					
					if(secondaryHeadloss < generalProperties.min_hl_perkm/1000 || secondaryHeadloss > generalProperties.max_hl_perkm/1000){
						linear = new Linear();
						linear.add(1, "l_"+i+"_"+fixedj+"_0");
						problem.add(new Constraint("", linear, "=", 0));
					}
					else
						atleastOneHeadlossCorrect = true;
					
					double primary_speed = Util.waterSpeed(flow, pipe.getDiameter());
					double secondary_speed = Util.waterSpeed(flow*secondaryFlowFactor, pipe.getDiameter());
					
					if(generalProperties.max_water_speed > 0 && primary_speed > generalProperties.max_water_speed){
						linear = new Linear();
						linear.add(1, "l_"+i+"_"+fixedj+"_1");
						problem.add(new Constraint("", linear, "=", 0));
					}
					else
						atleastOneSpeedCorrect = true;
					
					if(generalProperties.max_water_speed > 0 && secondary_speed > generalProperties.max_water_speed){
						linear = new Linear();
						linear.add(1, "l_"+i+"_"+fixedj+"_0");
						problem.add(new Constraint("", linear, "=", 0));
					}
					else
						atleastOneSpeedCorrect = true;
				}
			}
			else{
								
				int j=0;
				
				for(PipeCost entry : pipeCost){
					double flow = pipe.getFlow();
					double headloss = Util.HWheadLoss(flow, entry.getRoughness(), entry.getDiameter());
					double secondaryHeadloss = Util.HWheadLoss(flow * secondaryFlowFactor, entry.getRoughness(), entry.getDiameter());
					
					if(headloss < generalProperties.min_hl_perkm/1000 || headloss > generalProperties.max_hl_perkm/1000){
						linear = new Linear();
						linear.add(1, "l_"+i+"_"+j+"_1");
						problem.add(new Constraint("", linear, "=", 0));
					}
					else
						atleastOneHeadlossCorrect = true;
					
					if(secondaryHeadloss < generalProperties.min_hl_perkm/1000 || secondaryHeadloss > generalProperties.max_hl_perkm/1000){
						linear = new Linear();
						linear.add(1, "l_"+i+"_"+j+"_0");
						problem.add(new Constraint("", linear, "=", 0));
					}
					else
						atleastOneHeadlossCorrect = true;
					
					double primary_speed = Util.waterSpeed(flow, entry.getDiameter());
					double secondary_speed = Util.waterSpeed(flow*secondaryFlowFactor, entry.getDiameter());
					
					if(generalProperties.max_water_speed > 0 && primary_speed > generalProperties.max_water_speed){
						linear = new Linear();
						linear.add(1, "l_"+i+"_"+j+"_1");
						problem.add(new Constraint("", linear, "=", 0));
					}
					else
						atleastOneSpeedCorrect = true;
					
					if(generalProperties.max_water_speed > 0 && secondary_speed > generalProperties.max_water_speed){
						linear = new Linear();
						linear.add(1, "l_"+i+"_"+j+"_0");
						problem.add(new Constraint("", linear, "=", 0));
					}
					else
						atleastOneSpeedCorrect = true;
					
					j++;
				}
			}	
			if(!atleastOneHeadlossCorrect)
				throw new Exception("Headloss for pipe: " + pipe.getPipeID() + " cannot be set within the min/max headloss values.");
			if(!atleastOneSpeedCorrect)
				throw new Exception("Speed of water in pipe: " + pipe.getPipeID() + " cannot be set below the max value.");
		}
	}
		
	//set the headloss constraints in pipes
	private void setHeadLossConstraints() throws Exception{
		int j=0;
		for(Node node : nodes.values()){
			if(source!=node){
				Linear linear = new Linear();
				double valveloss = 0;
				for(Pipe pipe : node.getSourceToNodePipes()){
					int i = pipe.getPipeID();
					
					if(pipe.isAllowParallel()){
						j=0;
						for(PipeCost entry : pipeCost){
							// primary flow in case of parallel pipe
							double flow = pipe.getFlow() / (1 + (entry.getRoughness()/pipe.getRoughness())*Math.pow(entry.getDiameter()/pipe.getDiameter(), 4.87/1.852));
							double headloss = Util.HWheadLoss(pipe.getLength(), flow, pipe.getRoughness(), pipe.getDiameter());
							linear.add(headloss, "p_"+i+"_"+j);
							j++;
						}	
						// primary flow in case of no parallel pipe
						double flow = pipe.getFlow(); 
						double headloss = Util.HWheadLoss(pipe.getLength(), flow, pipe.getRoughness(), pipe.getDiameter());
						linear.add(headloss, "p_"+i);
					
					}
					else{
						if(pipe.getDiameter()!=0){
							j=0;
							for(PipeCost entry : pipeCost){
								double headloss = Util.HWheadLoss(pipe.getFlow(), pipe.getRoughness(), entry.getDiameter());
								linear.add(headloss, "l_"+i+"_"+j);
								j++;
							}
						}
						else{	
							j=0;
							for(PipeCost entry : pipeCost){
								double headloss = Util.HWheadLoss(pipe.getFlow(), entry.getRoughness(), entry.getDiameter());
								linear.add(headloss, "l_"+i+"_"+j);
								j++;
							}	
						}
					}
					
					//additional water head provided by pump 
					if(pumpGeneralProperties.pump_enabled && pipe.getAllowPump()){
						linear.add(-1, "pumphead_"+i);
					}
					
					//loss of water head due to pressure reducing valve
					valveloss += pipe.getValveSetting();
				}
				
				//System.out.println(linear + " " + (Node.getSourceHGL() - node.getElevation() - node.getResidualPressure()));
				//problem.add(linear, "<=", Node.getSourceHGL() - node.getElevation() - node.getResidualPressure());
				
//				//ESR height h_i for each node
//				if(node.getDemand()>0){
//					problem.setVarType("h_"+node.getNodeID(), Double.class);
//					problem.setVarLowerBound("h_"+node.getNodeID(), 0);
//					linear.add(1,"h_"+node.getNodeID());
//				}
				problem.add(new Constraint("", linear, "<=", source.getHead() - node.getElevation() - node.getResidualPressure() - valveloss));
			}
		}
	}

	//set headloss constraints in pipes
	//includes ESRs and thus primary and secondary network considerations
	private void setHeadLossConstraints_esr() throws Exception{
		Linear linear = new Linear();

		for(Node node : nodes.values()){
			int i = node.getNodeID();
			if(source==node){
				linear = new Linear();
				linear.add(1,"head_"+i);
				problem.add(new Constraint("", linear, "=", source.getHead()));
			}
			else{			
				// min_esr_height <= esr_i <= max_esr_height
				linear = new Linear();
				linear.add(minEsrHeight,"besr_"+i);
				linear.add(-1,"esr_"+i);
				problem.add(new Constraint("", linear, "<=",0));
								
				linear = new Linear();
				linear.add(maxEsrHeight,"besr_"+i);
				linear.add(-1,"esr_"+i);
				problem.add(new Constraint("", linear, ">=",0));
				
				linear = new Linear();
				linear.add(1,"besr_"+i);
				linear.add(-1,"s_"+i+"_"+i);
				problem.add(new Constraint("", linear, "=",0));
				
				linear = new Linear();
				linear.add(1,"head_"+i);
				linear.add(-1,"esr_"+i);
				problem.add(new Constraint("", linear, ">=", node.getElevation() + node.getResidualPressure()));
			}
		}		
		
		for(Pipe pipe : pipes.values()){
			linear = new Linear();
			int i = pipe.getPipeID();
			int j = 0;
			if(pipe.isAllowParallel()){
				for(PipeCost entry : pipeCost){
					// primary flow in case of parallel pipe
					double flow = pipe.getFlow() / (1 + (entry.getRoughness()/pipe.getRoughness())*Math.pow(entry.getDiameter()/pipe.getDiameter(), 4.87/1.852));
					double headloss = Util.HWheadLoss(pipe.getLength(), flow, pipe.getRoughness(), pipe.getDiameter());
					double headlossSecondary = Util.HWheadLoss(pipe.getLength(), flow * secondaryFlowFactor, pipe.getRoughness(), pipe.getDiameter());
					linear.add(headlossSecondary, "p_"+i+"_"+j);
					linear.add(headloss-headlossSecondary, "yp_"+i+"_"+j);
					j++;
				}	
				// primary flow in case of no parallel pipe
				double flow = pipe.getFlow(); 
				double headloss = Util.HWheadLoss(pipe.getLength(), flow, pipe.getRoughness(), pipe.getDiameter());
				double headlossSecondary = Util.HWheadLoss(pipe.getLength(), flow * secondaryFlowFactor, pipe.getRoughness(), pipe.getDiameter());
				linear.add(headlossSecondary, "p_"+i);
				linear.add(headloss-headlossSecondary, "yp_"+i);					
			}
			else{
				if(pipe.getDiameter()!=0){
					j=0;
					for(PipeCost entry : pipeCost){
						double headloss = Util.HWheadLoss(pipe.getFlow(), pipe.getRoughness(), entry.getDiameter());
						double headlossSecondary = Util.HWheadLoss(pipe.getFlow() * secondaryFlowFactor, pipe.getRoughness(), entry.getDiameter());
						linear.add(headlossSecondary, "l_"+i+"_"+j);
						linear.add(headloss-headlossSecondary, "y_"+i+"_"+j);
						j++;
					}
				}
				else{	
					j=0;
					for(PipeCost entry : pipeCost){
						double headloss = Util.HWheadLoss(pipe.getFlow(), entry.getRoughness(), entry.getDiameter());
						double headlossSecondary = Util.HWheadLoss(pipe.getFlow() * secondaryFlowFactor, entry.getRoughness(), entry.getDiameter());
						linear.add(headlossSecondary, "l_"+i+"_"+j);
						linear.add(headloss-headlossSecondary, "y_"+i+"_"+j);
						j++;
					}	
				}
			}
			Node startNode = pipe.getStartNode();
			int startid = startNode.getNodeID();
			int destinationid = pipe.getEndNode().getNodeID();
			
			linear.add(-1,"headloss_"+i);
			problem.add(new Constraint("", linear, "=",0));
			
			
			// head_source = head_source if demand_source = 0
			// head_source = head_elevation + esr_height if s_source_source = 1 and f_pipe = 0
			// else head_source = head_source
			// introduce head_i_j : is the start node, j is the id of the pipe
			// k_i_j = s_i_i && !f_j
			// head_i_j = k_i_j * (elevation_i + esr_i - head_i) + head_i
			// ehead_i = elevation_i + esr_i - head_i
			// yhead_i_j = k_i_j * ehead_i
			
			if(startNode.getDemand()>0){				
				// head_i_j = yhead_i_j + head_i
				linear = new Linear();
				linear.add(1,"yhead_"+startid+"_"+i);
				linear.add(1,"head_"+startid);
				linear.add(-1,"head_"+startid+"_"+i);
				problem.add(new Constraint("", linear, "=",0));
				
				// ehead_i = elevation_i + esr_i - head_i
				linear = new Linear();
				linear.add(1,"ehead_"+startid);
				linear.add(1,"head_"+startid);
				linear.add(-1,"esr_"+startid);
				problem.add(new Constraint("", linear, "=",startNode.getElevation()));
				
				//to capture yhead_startid_i = k_startid_i * ehead_startid
				// ehead range assumed to be -10000,10000
				int minEHead = -10000;
				int maxEHead = 10000;
						
				linear = new Linear();
				linear.add(1,"yhead_"+startid+"_"+i);
				problem.add(new Constraint("", linear, "<=",maxEHead));
				
				linear = new Linear();
				linear.add(1,"yhead_"+startid+"_"+i);
				problem.add(new Constraint("", linear, ">=",minEHead));
				
				linear = new Linear();
				linear.add(-1,"yhead_"+startid+"_"+i);
				linear.add(maxEHead,"k_"+startid+"_"+i);
				problem.add(new Constraint("", linear, ">=",0));
				
				linear = new Linear();
				linear.add(-1,"yhead_"+startid+"_"+i);
				linear.add(minEHead,"k_"+startid+"_"+i);
				problem.add(new Constraint("", linear, "<=",0));
				
				linear = new Linear();
				linear.add(1,"ehead_"+startid);
				linear.add(maxEHead,"k_"+startid+"_"+i);
				linear.add(-1,"yhead_"+startid+"_"+i);
				problem.add(new Constraint("", linear, "<=",maxEHead));
				
				linear = new Linear();
				linear.add(1,"ehead_"+startid);
				linear.add(minEHead,"k_"+startid+"_"+i);
				linear.add(-1,"yhead_"+startid+"_"+i);
				problem.add(new Constraint("", linear, ">=",minEHead));
				
				linear = new Linear();
				linear.add(-1,"ehead_"+startid);
				linear.add(maxEHead,"k_"+startid+"_"+i);
				linear.add(1,"yhead_"+startid+"_"+i);
				problem.add(new Constraint("", linear, "<=",maxEHead));
				
				
				//to capture k_startid_i = s_startid_startid ^ !f_i
				
				linear = new Linear();
				linear.add(1,"s_"+startid+"_"+startid);
				linear.add(-1,"k_"+startid+"_"+i);
				problem.add(new Constraint("", linear, ">=",0));
				
				linear = new Linear();
				linear.add(1,"f_"+i);
				linear.add(1,"k_"+startid+"_"+i);
				problem.add(new Constraint("", linear, "<=",1));
				
				linear = new Linear();
				linear.add(1,"f_"+i);
				linear.add(1,"k_"+startid+"_"+i);
				linear.add(-1,"s_"+startid+"_"+startid);
				problem.add(new Constraint("", linear, ">=",0));
				
				//head_destination = head_i_j - headloss
				linear = new Linear();
				linear.add(1,"headloss_"+i);
				linear.add(1,"head_"+destinationid);
				linear.add(-1,"head_"+startid+"_"+i);
				problem.add(new Constraint("", linear, "=",0));
			}
			else{
				// head_destination = head_source - headloss
				linear = new Linear();
				linear.add(1,"headloss_"+i);
				linear.add(1,"head_"+destinationid);
				linear.add(-1,"head_"+startid);
				problem.add(new Constraint("", linear, "=",0));
			}
		}
	}
	
	//set headloss constraints in pipes
	//includes ESRs and thus primary and secondary network considerations
	//also allows nodes with zero demands to have ESRs
	//allows removal of certain nodes from ESR consideration
	private void setHeadLossConstraints_esr_gen2() throws Exception{
		Linear linear = new Linear();

		for(Node node : nodes.values()){
			int i = node.getNodeID();
			if(source==node){
				linear = new Linear();
				linear.add(1,"head_"+i);
				problem.add(new Constraint("", linear, "=", source.getHead()));
			}
			else{
				linear = new Linear();
				linear.add(1,"head_"+i);
				linear.add(-1,"esr_"+i);
				problem.add(new Constraint("", linear, ">=", node.getElevation() + node.getResidualPressure()));
			
			}
				
			if(!node.getAllowESR()){
				linear = new Linear();
				linear.add(1,"d_"+i);
				problem.add(new Constraint("", linear, "=", 0));
			
			}
			
			// min_esr_height <= esr_i <= max_esr_height
			linear = new Linear();
			linear.add(minEsrHeight,"besr_"+i);
			linear.add(-1,"esr_"+i);
			problem.add(new Constraint("", linear, "<=",0));
							
			linear = new Linear();
			linear.add(maxEsrHeight,"besr_"+i);
			linear.add(-1,"esr_"+i);
			problem.add(new Constraint("", linear, ">=",0));
			
			linear = new Linear();
			linear.add(1,"besr_"+i);
			linear.add(-1,"s_"+i+"_"+i);
			problem.add(new Constraint("", linear, "=",0));
			
		}		
		
		for(Pipe pipe : pipes.values()){
			linear = new Linear();
			int i = pipe.getPipeID();
			int j = 0;
			if(pipe.isAllowParallel()){
				for(PipeCost entry : pipeCost){
					// primary flow in case of parallel pipe
					double flow = pipe.getFlow() / (1 + (entry.getRoughness()/pipe.getRoughness())*Math.pow(entry.getDiameter()/pipe.getDiameter(), 4.87/1.852));
					double headloss = Util.HWheadLoss(pipe.getLength(), flow, pipe.getRoughness(), pipe.getDiameter());
					double headlossSecondary = Util.HWheadLoss(pipe.getLength(), flow * secondaryFlowFactor, pipe.getRoughness(), pipe.getDiameter());
					linear.add(headlossSecondary, "p_"+i+"_"+j);
					linear.add(headloss-headlossSecondary, "yp_"+i+"_"+j);
					j++;
				}	
				// primary flow in case of no parallel pipe
				double flow = pipe.getFlow(); 
				double headloss = Util.HWheadLoss(pipe.getLength(), flow, pipe.getRoughness(), pipe.getDiameter());
				double headlossSecondary = Util.HWheadLoss(pipe.getLength(), flow * secondaryFlowFactor, pipe.getRoughness(), pipe.getDiameter());
				linear.add(headlossSecondary, "p_"+i);
				linear.add(headloss-headlossSecondary, "yp_"+i);					
			}
			else{
				if(pipe.getDiameter()!=0){
					j=0;
					for(PipeCost entry : pipeCost){
						double headloss = Util.HWheadLoss(pipe.getFlow(), pipe.getRoughness(), entry.getDiameter());
						double headlossSecondary = Util.HWheadLoss(pipe.getFlow() * secondaryFlowFactor, pipe.getRoughness(), entry.getDiameter());
						linear.add(headlossSecondary, "l_"+i+"_"+j);
						linear.add(headloss-headlossSecondary, "y_"+i+"_"+j);
						j++;
					}
				}
				else{	
					j=0;
					for(PipeCost entry : pipeCost){
						double headloss = Util.HWheadLoss(pipe.getFlow(), entry.getRoughness(), entry.getDiameter());
						double headlossSecondary = Util.HWheadLoss(pipe.getFlow() * secondaryFlowFactor, entry.getRoughness(), entry.getDiameter());
						linear.add(headlossSecondary, "l_"+i+"_"+j);
						linear.add(headloss-headlossSecondary, "y_"+i+"_"+j);
						j++;
					}	
				}
			}
			Node startNode = pipe.getStartNode();
			int startid = startNode.getNodeID();
			int destinationid = pipe.getEndNode().getNodeID();
			
			linear.add(-1,"headloss_"+i);
			problem.add(new Constraint("", linear, "=",0));
			
			
			// head_source = head_source if demand_source = 0
			// head_source = head_elevation + esr_height if s_source_source = 1 and f_pipe = 0
			// else head_source = head_source
			// introduce head_i_j : is the start node, j is the id of the pipe
			// k_i_j = s_i_i && !f_j
			// head_i_j = k_i_j * (elevation_i + esr_i - head_i) + head_i
			// ehead_i = elevation_i + esr_i - head_i
			// yhead_i_j = k_i_j * ehead_i
					
			if(startNode.getAllowESR()){
				// head_i_j = yhead_i_j + head_i
				linear = new Linear();
				linear.add(1,"yhead_"+startid+"_"+i);
				linear.add(1,"head_"+startid);
				linear.add(-1,"head_"+startid+"_"+i);
				problem.add(new Constraint("", linear, "=",0));
				
				// ehead_i = elevation_i + esr_i - head_i
				linear = new Linear();
				linear.add(1,"ehead_"+startid);
				linear.add(1,"head_"+startid);
				linear.add(-1,"esr_"+startid);
				problem.add(new Constraint("", linear, "=",startNode.getElevation()));
				
				//to capture yhead_startid_i = k_startid_i * ehead_startid
				// ehead range assumed to be -10000,10000
				int minEHead = -10000;
				int maxEHead = 10000;
						
				linear = new Linear();
				linear.add(1,"yhead_"+startid+"_"+i);
				problem.add(new Constraint("", linear, "<=",maxEHead));
				
				linear = new Linear();
				linear.add(1,"yhead_"+startid+"_"+i);
				problem.add(new Constraint("", linear, ">=",minEHead));
				
				linear = new Linear();
				linear.add(-1,"yhead_"+startid+"_"+i);
				linear.add(maxEHead,"k_"+startid+"_"+i);
				problem.add(new Constraint("", linear, ">=",0));
				
				linear = new Linear();
				linear.add(-1,"yhead_"+startid+"_"+i);
				linear.add(minEHead,"k_"+startid+"_"+i);
				problem.add(new Constraint("", linear, "<=",0));
				
				linear = new Linear();
				linear.add(1,"ehead_"+startid);
				linear.add(maxEHead,"k_"+startid+"_"+i);
				linear.add(-1,"yhead_"+startid+"_"+i);
				problem.add(new Constraint("", linear, "<=",maxEHead));
				
				linear = new Linear();
				linear.add(1,"ehead_"+startid);
				linear.add(minEHead,"k_"+startid+"_"+i);
				linear.add(-1,"yhead_"+startid+"_"+i);
				problem.add(new Constraint("", linear, ">=",minEHead));
				
				linear = new Linear();
				linear.add(-1,"ehead_"+startid);
				linear.add(maxEHead,"k_"+startid+"_"+i);
				linear.add(1,"yhead_"+startid+"_"+i);
				problem.add(new Constraint("", linear, "<=",maxEHead));
				
				
				//to capture k_startid_i = s_startid_startid ^ !f_i
				
				linear = new Linear();
				linear.add(1,"s_"+startid+"_"+startid);
				linear.add(-1,"k_"+startid+"_"+i);
				problem.add(new Constraint("", linear, ">=",0));
				
				linear = new Linear();
				linear.add(1,"f_"+i);
				linear.add(1,"k_"+startid+"_"+i);
				problem.add(new Constraint("", linear, "<=",1));
				
				linear = new Linear();
				linear.add(1,"f_"+i);
				linear.add(1,"k_"+startid+"_"+i);
				linear.add(-1,"s_"+startid+"_"+startid);
				problem.add(new Constraint("", linear, ">=",0));
				
				//head_destination = head_i_j - headloss
				linear = new Linear();
				linear.add(1,"headloss_"+i);
				linear.add(1,"head_"+destinationid);
				linear.add(-1,"head_"+startid+"_"+i);
				problem.add(new Constraint("", linear, "=",0));
			}
			else{
				// head_destination = head_source - headloss
				linear = new Linear();
				linear.add(1,"headloss_"+i);
				linear.add(1,"head_"+destinationid);
				linear.add(-1,"head_"+startid);
				problem.add(new Constraint("", linear, "=",0));
			}
		}
	}
	
	//set headloss constraints in pipes
	//includes ESRs and thus primary and secondary network considerations
	//also allows nodes with zero demands to have ESRs
	//updated to use l_i_j_k instead of l_i_j
	private void setHeadLossConstraints_esr_gen3() throws Exception{
		Linear linear = new Linear();

		for(Node node : nodes.values()){
			int i = node.getNodeID();
			if(source==node){
				linear = new Linear();
				linear.add(1,"head_"+i);
				problem.add(new Constraint("", linear, "=", source.getHead()));
			}
			else{
				linear = new Linear();
				linear.add(1,"head_"+i);
				linear.add(-1,"esr_"+i);
				problem.add(new Constraint("", linear, ">=", node.getElevation() + node.getResidualPressure()));
			
			}
				
			if(!node.getAllowESR()){
				linear = new Linear();
				linear.add(1,"d_"+i);
				problem.add(new Constraint("", linear, "=", 0));
			
			}
			
			// min_esr_height <= esr_i <= max_esr_height
			linear = new Linear();
			linear.add(minEsrHeight,"besr_"+i);
			linear.add(-1,"esr_"+i);
			problem.add(new Constraint("", linear, "<=",0));
							
			linear = new Linear();
			linear.add(maxEsrHeight,"besr_"+i);
			linear.add(-1,"esr_"+i);
			problem.add(new Constraint("", linear, ">=",0));
			
			linear = new Linear();
			linear.add(1,"besr_"+i);
			linear.add(-1,"s_"+i+"_"+i);
			problem.add(new Constraint("", linear, "=",0));
			
		}		
		
		for(Pipe pipe : pipes.values()){
			linear = new Linear();
			int i = pipe.getPipeID();
			int j = 0;
			if(pipe.isAllowParallel()){
				for(PipeCost entry : pipeCost){
					// primary flow in case of parallel pipe
					double flow = pipe.getFlow() / (1 + (entry.getRoughness()/pipe.getRoughness())*Math.pow(entry.getDiameter()/pipe.getDiameter(), 4.87/1.852));
					double headloss = Util.HWheadLoss(pipe.getLength(), flow, pipe.getRoughness(), pipe.getDiameter());
					double headlossSecondary = Util.HWheadLoss(pipe.getLength(), flow * secondaryFlowFactor, pipe.getRoughness(), pipe.getDiameter());
					linear.add(headlossSecondary, "p_"+i+"_"+j+"_0");
					linear.add(headloss, "p_"+i+"_"+j+"_1");
					j++;
				}	
				// primary flow in case of no parallel pipe
				double flow = pipe.getFlow(); 
				double headloss = Util.HWheadLoss(pipe.getLength(), flow, pipe.getRoughness(), pipe.getDiameter());
				double headlossSecondary = Util.HWheadLoss(pipe.getLength(), flow * secondaryFlowFactor, pipe.getRoughness(), pipe.getDiameter());
				linear.add(headlossSecondary, "p_"+i+"_0");
				linear.add(headloss, "p_"+i+"_1");					
			}
			else{
				if(pipe.getDiameter()!=0){
					j=0;
					for(PipeCost entry : pipeCost){
						double headloss = Util.HWheadLoss(pipe.getFlow(), pipe.getRoughness(), entry.getDiameter());
						double headlossSecondary = Util.HWheadLoss(pipe.getFlow() * secondaryFlowFactor, pipe.getRoughness(), entry.getDiameter());
						linear.add(headlossSecondary, "l_"+i+"_"+j+"_0");
						linear.add(headloss, "l_"+i+"_"+j+"_1");
						j++;
					}
				}
				else{	
					j=0;
					for(PipeCost entry : pipeCost){
						double headloss = Util.HWheadLoss(pipe.getFlow(), entry.getRoughness(), entry.getDiameter());
						double headlossSecondary = Util.HWheadLoss(pipe.getFlow() * secondaryFlowFactor, entry.getRoughness(), entry.getDiameter());
						linear.add(headlossSecondary, "l_"+i+"_"+j+"_0");
						linear.add(headloss, "l_"+i+"_"+j+"_1");
						j++;
					}	
				}
			}
			linear.add(-1,"headloss_"+i);
			problem.add(new Constraint("", linear, "=",0));
			
			
			Node startNode = pipe.getStartNode();
			int startid = startNode.getNodeID();
			int destinationid = pipe.getEndNode().getNodeID();
						
			// head_source = head_source if demand_source = 0
			// head_source = head_elevation + esr_height if s_source_source = 1 and f_pipe = 0
			// else head_source = head_source
			// introduce head_i_j : is the start node, j is the id of the pipe
			// k_i_j = s_i_i && !f_j
			// head_i_j = k_i_j * (elevation_i + esr_i - head_i) + head_i
			// ehead_i = elevation_i + esr_i - head_i
			// yhead_i_j = k_i_j * ehead_i
					
			if(startNode.getAllowESR()){
				// head_i_j = yhead_i_j + head_i
				linear = new Linear();
				linear.add(1,"yhead_"+startid+"_"+i);
				linear.add(1,"head_"+startid);
				linear.add(-1,"head_"+startid+"_"+i);
				problem.add(new Constraint("", linear, "=",0));
				
				// ehead_i = elevation_i + esr_i - head_i
				linear = new Linear();
				linear.add(1,"ehead_"+startid);
				linear.add(1,"head_"+startid);
				linear.add(-1,"esr_"+startid);
				problem.add(new Constraint("", linear, "=",startNode.getElevation()));
				
				//to capture yhead_startid_i = k_startid_i * ehead_startid
				// ehead range assumed to be -10000,10000
				int minEHead = -10000;
				int maxEHead = 10000;
						
				linear = new Linear();
				linear.add(1,"yhead_"+startid+"_"+i);
				problem.add(new Constraint("", linear, "<=",maxEHead));
				
				linear = new Linear();
				linear.add(1,"yhead_"+startid+"_"+i);
				problem.add(new Constraint("", linear, ">=",minEHead));
				
				linear = new Linear();
				linear.add(-1,"yhead_"+startid+"_"+i);
				linear.add(maxEHead,"k_"+startid+"_"+i);
				problem.add(new Constraint("", linear, ">=",0));
				
				linear = new Linear();
				linear.add(-1,"yhead_"+startid+"_"+i);
				linear.add(minEHead,"k_"+startid+"_"+i);
				problem.add(new Constraint("", linear, "<=",0));
				
				linear = new Linear();
				linear.add(1,"ehead_"+startid);
				linear.add(maxEHead,"k_"+startid+"_"+i);
				linear.add(-1,"yhead_"+startid+"_"+i);
				problem.add(new Constraint("", linear, "<=",maxEHead));
				
				linear = new Linear();
				linear.add(1,"ehead_"+startid);
				linear.add(minEHead,"k_"+startid+"_"+i);
				linear.add(-1,"yhead_"+startid+"_"+i);
				problem.add(new Constraint("", linear, ">=",minEHead));
				
				linear = new Linear();
				linear.add(-1,"ehead_"+startid);
				linear.add(maxEHead,"k_"+startid+"_"+i);
				linear.add(1,"yhead_"+startid+"_"+i);
				problem.add(new Constraint("", linear, "<=",maxEHead));
				
				
				//to capture k_startid_i = s_startid_startid ^ !f_i
				
				linear = new Linear();
				linear.add(1,"s_"+startid+"_"+startid);
				linear.add(-1,"k_"+startid+"_"+i);
				problem.add(new Constraint("", linear, ">=",0));
				
				linear = new Linear();
				linear.add(1,"f_"+i);
				linear.add(1,"k_"+startid+"_"+i);
				problem.add(new Constraint("", linear, "<=",1));
				
				linear = new Linear();
				linear.add(1,"f_"+i);
				linear.add(1,"k_"+startid+"_"+i);
				linear.add(-1,"s_"+startid+"_"+startid);
				problem.add(new Constraint("", linear, ">=",0));
				
				//head_destination = head_i_j - headloss
				linear = new Linear();
				linear.add(1,"headloss_"+i);
				linear.add(1,"head_"+destinationid);
				linear.add(-1,"head_"+startid+"_"+i);
				problem.add(new Constraint("", linear, "=",0));
			}
			else{
				// head_destination = head_source - headloss
				linear = new Linear();
				linear.add(1,"headloss_"+i);
				linear.add(1,"head_"+destinationid);
				linear.add(-1,"head_"+startid);
				problem.add(new Constraint("", linear, "=",0));
			}
		}
	}
	
	//set headloss constraints in pipes
	//includes ESRs and thus primary and secondary network considerations
	//also allows nodes with zero demands to have ESRs
	//updated to use l_i_j_k instead of l_i_j
	//includes pumps
	private void setHeadLossConstraints_esr_gen4() throws Exception{
		Linear linear = new Linear();

		for(Node node : nodes.values()){
			int i = node.getNodeID();
			if(source==node){
				linear = new Linear();
				linear.add(1,"head_"+i);
				problem.add(new Constraint("", linear, "=", source.getHead()));
			}
			else{
				linear = new Linear();
				linear.add(1,"head_"+i);
				linear.add(-1,"esr_"+i);
				problem.add(new Constraint("", linear, ">=", node.getElevation() + node.getResidualPressure()));
			
			}
				
			if(!node.getAllowESR()){
				linear = new Linear();
				linear.add(1,"d_"+i);
				problem.add(new Constraint("", linear, "=", 0));
			
			}
			
			// min_esr_height <= esr_i <= max_esr_height
			linear = new Linear();
			linear.add(minEsrHeight,"besr_"+i);
			linear.add(-1,"esr_"+i);
			problem.add(new Constraint("", linear, "<=",0));
							
			linear = new Linear();
			linear.add(maxEsrHeight,"besr_"+i);
			linear.add(-1,"esr_"+i);
			problem.add(new Constraint("", linear, ">=",0));
			
			linear = new Linear();
			linear.add(1,"besr_"+i);
			linear.add(-1,"s_"+i+"_"+i);
			problem.add(new Constraint("", linear, "=",0));
			
		}		
		
		for(Pipe pipe : pipes.values()){
			linear = new Linear();
			int i = pipe.getPipeID();
			int j = 0;
			if(pipe.isAllowParallel()){
				for(PipeCost entry : pipeCost){
					// primary flow in case of parallel pipe
					double flow = pipe.getFlow() / (1 + (entry.getRoughness()/pipe.getRoughness())*Math.pow(entry.getDiameter()/pipe.getDiameter(), 4.87/1.852));
					double headloss = Util.HWheadLoss(pipe.getLength(), flow, pipe.getRoughness(), pipe.getDiameter());
					double headlossSecondary = Util.HWheadLoss(pipe.getLength(), flow * secondaryFlowFactor, pipe.getRoughness(), pipe.getDiameter());
					linear.add(headlossSecondary, "p_"+i+"_"+j+"_0");
					linear.add(headloss, "p_"+i+"_"+j+"_1");
					j++;
				}	
				// primary flow in case of no parallel pipe
				double flow = pipe.getFlow(); 
				double headloss = Util.HWheadLoss(pipe.getLength(), flow, pipe.getRoughness(), pipe.getDiameter());
				double headlossSecondary = Util.HWheadLoss(pipe.getLength(), flow * secondaryFlowFactor, pipe.getRoughness(), pipe.getDiameter());
				linear.add(headlossSecondary, "p_"+i+"_0");
				linear.add(headloss, "p_"+i+"_1");					
			}
			else{
				if(pipe.getDiameter()!=0){
					j=0;
					for(PipeCost entry : pipeCost){
						double headloss = Util.HWheadLoss(pipe.getFlow(), pipe.getRoughness(), entry.getDiameter());
						double headlossSecondary = Util.HWheadLoss(pipe.getFlow() * secondaryFlowFactor, pipe.getRoughness(), entry.getDiameter());
						linear.add(headlossSecondary, "l_"+i+"_"+j+"_0");
						linear.add(headloss, "l_"+i+"_"+j+"_1");
						j++;
					}
				}
				else{	
					j=0;
					for(PipeCost entry : pipeCost){
						double headloss = Util.HWheadLoss(pipe.getFlow(), entry.getRoughness(), entry.getDiameter());
						double headlossSecondary = Util.HWheadLoss(pipe.getFlow() * secondaryFlowFactor, entry.getRoughness(), entry.getDiameter());
						linear.add(headlossSecondary, "l_"+i+"_"+j+"_0");
						linear.add(headloss, "l_"+i+"_"+j+"_1");
						j++;
					}	
				}
			}
			
			if(pumpGeneralProperties.pump_enabled && pipe.getAllowPump()){
				linear.add(-1, "pumphead_"+i);
			}
			
			linear.add(-1,"headloss_"+i);
			problem.add(new Constraint("", linear, "=",-1*pipe.getValveSetting()));
					
			Node startNode = pipe.getStartNode();
			int startid = startNode.getNodeID();
			int destinationid = pipe.getEndNode().getNodeID();
						
			// head_source = head_source if demand_source = 0
			// head_source = head_elevation + esr_height if s_source_source = 1 and f_pipe = 0
			// else head_source = head_source
			// introduce head_i_j : is the start node, j is the id of the pipe
			// k_i_j = s_i_i && !f_j
			// head_i_j = k_i_j * (elevation_i + esr_i - head_i) + head_i
			// ehead_i = elevation_i + esr_i - head_i
			// yhead_i_j = k_i_j * ehead_i
					
			if(startNode.getAllowESR()){
				// head_i_j = yhead_i_j + head_i
				linear = new Linear();
				linear.add(1,"yhead_"+startid+"_"+i);
				linear.add(1,"head_"+startid);
				linear.add(-1,"head_"+startid+"_"+i);
				problem.add(new Constraint("", linear, "=",0));
				
				// ehead_i = elevation_i + esr_i - head_i
				linear = new Linear();
				linear.add(1,"ehead_"+startid);
				linear.add(1,"head_"+startid);
				linear.add(-1,"esr_"+startid);
				problem.add(new Constraint("", linear, "=",startNode.getElevation()));
				
				//to capture yhead_startid_i = k_startid_i * ehead_startid
				// ehead range assumed to be -10000,10000
				int minEHead = -10000;
				int maxEHead = 10000;
						
				linear = new Linear();
				linear.add(1,"yhead_"+startid+"_"+i);
				problem.add(new Constraint("", linear, "<=",maxEHead));
				
				linear = new Linear();
				linear.add(1,"yhead_"+startid+"_"+i);
				problem.add(new Constraint("", linear, ">=",minEHead));
				
				linear = new Linear();
				linear.add(-1,"yhead_"+startid+"_"+i);
				linear.add(maxEHead,"k_"+startid+"_"+i);
				problem.add(new Constraint("", linear, ">=",0));
				
				linear = new Linear();
				linear.add(-1,"yhead_"+startid+"_"+i);
				linear.add(minEHead,"k_"+startid+"_"+i);
				problem.add(new Constraint("", linear, "<=",0));
				
				linear = new Linear();
				linear.add(1,"ehead_"+startid);
				linear.add(maxEHead,"k_"+startid+"_"+i);
				linear.add(-1,"yhead_"+startid+"_"+i);
				problem.add(new Constraint("", linear, "<=",maxEHead));
				
				linear = new Linear();
				linear.add(1,"ehead_"+startid);
				linear.add(minEHead,"k_"+startid+"_"+i);
				linear.add(-1,"yhead_"+startid+"_"+i);
				problem.add(new Constraint("", linear, ">=",minEHead));
				
				linear = new Linear();
				linear.add(-1,"ehead_"+startid);
				linear.add(maxEHead,"k_"+startid+"_"+i);
				linear.add(1,"yhead_"+startid+"_"+i);
				problem.add(new Constraint("", linear, "<=",maxEHead));
				
				
				//to capture k_startid_i = s_startid_startid ^ !f_i
				
				linear = new Linear();
				linear.add(1,"s_"+startid+"_"+startid);
				linear.add(-1,"k_"+startid+"_"+i);
				problem.add(new Constraint("", linear, ">=",0));
				
				linear = new Linear();
				linear.add(1,"f_"+i);
				linear.add(1,"k_"+startid+"_"+i);
				problem.add(new Constraint("", linear, "<=",1));
				
				linear = new Linear();
				linear.add(1,"f_"+i);
				linear.add(1,"k_"+startid+"_"+i);
				linear.add(-1,"s_"+startid+"_"+startid);
				problem.add(new Constraint("", linear, ">=",0));
				
				//head_destination = head_i_j - headloss
				linear = new Linear();
				linear.add(1,"headloss_"+i);
				linear.add(1,"head_"+destinationid);
				linear.add(-1,"head_"+startid+"_"+i);
				problem.add(new Constraint("", linear, "=",0));
			}
			else{
				// head_destination = head_source - headloss
				linear = new Linear();
				linear.add(1,"headloss_"+i);
				linear.add(1,"head_"+destinationid);
				linear.add(-1,"head_"+startid);
				problem.add(new Constraint("", linear, "=",0));
			}
		}
	}
	
	//set headloss constraints in pipes
	//removed usage of s_i_j variables
	private void setHeadLossConstraints_esr_gen5() throws Exception{
		Linear linear = new Linear();

		for(Node node : nodes.values()){
			int i = node.getNodeID();
			if(source==node){
				linear = new Linear();
				linear.add(1,"head_"+i);
				problem.add(new Constraint("", linear, "=", source.getHead()));
			}
			else{
				linear = new Linear();
				linear.add(1,"head_"+i);
				linear.add(-1,"esr_"+i);
				problem.add(new Constraint("", linear, ">=", node.getElevation() + node.getResidualPressure()));
			
			}
				
			if(!node.getAllowESR()){
				linear = new Linear();
				linear.add(1,"d_"+i);
				problem.add(new Constraint("", linear, "=", 0));
			
			}
			
			// min_esr_height <= esr_i <= max_esr_height
			linear = new Linear();
			linear.add(minEsrHeight,"besr_"+i);
			linear.add(-1,"esr_"+i);
			problem.add(new Constraint("", linear, "<=",0));
							
			linear = new Linear();
			linear.add(maxEsrHeight,"besr_"+i);
			linear.add(-1,"esr_"+i);
			problem.add(new Constraint("", linear, ">=",0));			
		}		
		
		for(Pipe pipe : pipes.values()){
			linear = new Linear();
			int i = pipe.getPipeID();
			int j = 0;
			if(pipe.isAllowParallel()){
				for(PipeCost entry : pipeCost){
					// primary flow in case of parallel pipe
					double flow = pipe.getFlow() / (1 + (entry.getRoughness()/pipe.getRoughness())*Math.pow(entry.getDiameter()/pipe.getDiameter(), 4.87/1.852));
					double headloss = Util.HWheadLoss(pipe.getLength(), flow, pipe.getRoughness(), pipe.getDiameter());
					double headlossSecondary = Util.HWheadLoss(pipe.getLength(), flow * secondaryFlowFactor, pipe.getRoughness(), pipe.getDiameter());
					linear.add(headlossSecondary, "p_"+i+"_"+j+"_0");
					linear.add(headloss, "p_"+i+"_"+j+"_1");
					j++;
				}	
				// primary flow in case of no parallel pipe
				double flow = pipe.getFlow(); 
				double headloss = Util.HWheadLoss(pipe.getLength(), flow, pipe.getRoughness(), pipe.getDiameter());
				double headlossSecondary = Util.HWheadLoss(pipe.getLength(), flow * secondaryFlowFactor, pipe.getRoughness(), pipe.getDiameter());
				linear.add(headlossSecondary, "p_"+i+"_0");
				linear.add(headloss, "p_"+i+"_1");					
			}
			else{
				if(pipe.getDiameter()!=0){
					j=0;
					for(PipeCost entry : pipeCost){
						double headloss = Util.HWheadLoss(pipe.getFlow(), pipe.getRoughness(), entry.getDiameter());
						double headlossSecondary = Util.HWheadLoss(pipe.getFlow() * secondaryFlowFactor, pipe.getRoughness(), entry.getDiameter());
						linear.add(headlossSecondary, "l_"+i+"_"+j+"_0");
						linear.add(headloss, "l_"+i+"_"+j+"_1");
						j++;
					}
				}
				else{	
					j=0;
					for(PipeCost entry : pipeCost){
						double headloss = Util.HWheadLoss(pipe.getFlow(), entry.getRoughness(), entry.getDiameter());
						double headlossSecondary = Util.HWheadLoss(pipe.getFlow() * secondaryFlowFactor, entry.getRoughness(), entry.getDiameter());
						linear.add(headlossSecondary, "l_"+i+"_"+j+"_0");
						linear.add(headloss, "l_"+i+"_"+j+"_1");
						j++;
					}	
				}
			}
			
			if(pumpGeneralProperties.pump_enabled && pipe.getAllowPump()){
				linear.add(-1, "pumphead_"+i);
			}
			
			linear.add(-1,"headloss_"+i);
			problem.add(new Constraint("", linear, "=",-1*pipe.getValveSetting()));
					
			Node startNode = pipe.getStartNode();
			int startid = startNode.getNodeID();
			int destinationid = pipe.getEndNode().getNodeID();
						
			// head_source = head_source if demand_source = 0
			// head_source = head_elevation + esr_height if s_source_source = 1 and f_pipe = 0
			// else head_source = head_source
			// introduce head_i_j : i is the start node, j is the id of the pipe
			// k_i_j = besr_i && !f_j
			// head_i_j = k_i_j * (elevation_i + esr_i - head_i) + head_i
			// ehead_i = elevation_i + esr_i - head_i
			// yhead_i_j = k_i_j * ehead_i
					
			if(startNode.getAllowESR()){
				// head_i_j = yhead_i_j + head_i
				linear = new Linear();
				linear.add(1,"yhead_"+startid+"_"+i);
				linear.add(1,"head_"+startid);
				linear.add(-1,"head_"+startid+"_"+i);
				problem.add(new Constraint("", linear, "=",0));
				
				// ehead_i = elevation_i + esr_i - head_i
				linear = new Linear();
				linear.add(1,"ehead_"+startid);
				linear.add(1,"head_"+startid);
				linear.add(-1,"esr_"+startid);
				problem.add(new Constraint("", linear, "=",startNode.getElevation()));
				
				//to capture yhead_startid_i = k_startid_i * ehead_startid
				// ehead range assumed to be -10000,10000
				int minEHead = -10000;
				int maxEHead = 10000;
						
				linear = new Linear();
				linear.add(1,"yhead_"+startid+"_"+i);
				problem.add(new Constraint("", linear, "<=",maxEHead));
				
				linear = new Linear();
				linear.add(1,"yhead_"+startid+"_"+i);
				problem.add(new Constraint("", linear, ">=",minEHead));
				
				linear = new Linear();
				linear.add(-1,"yhead_"+startid+"_"+i);
				linear.add(maxEHead,"k_"+startid+"_"+i);
				problem.add(new Constraint("", linear, ">=",0));
				
				linear = new Linear();
				linear.add(-1,"yhead_"+startid+"_"+i);
				linear.add(minEHead,"k_"+startid+"_"+i);
				problem.add(new Constraint("", linear, "<=",0));
				
				linear = new Linear();
				linear.add(1,"ehead_"+startid);
				linear.add(maxEHead,"k_"+startid+"_"+i);
				linear.add(-1,"yhead_"+startid+"_"+i);
				problem.add(new Constraint("", linear, "<=",maxEHead));
				
				linear = new Linear();
				linear.add(1,"ehead_"+startid);
				linear.add(minEHead,"k_"+startid+"_"+i);
				linear.add(-1,"yhead_"+startid+"_"+i);
				problem.add(new Constraint("", linear, ">=",minEHead));
				
				linear = new Linear();
				linear.add(-1,"ehead_"+startid);
				linear.add(maxEHead,"k_"+startid+"_"+i);
				linear.add(1,"yhead_"+startid+"_"+i);
				problem.add(new Constraint("", linear, "<=",maxEHead));
				
				
				//to capture k_startid_i = besr_startid ^ !f_i
				
				linear = new Linear();
				linear.add(1,"besr_"+startid);
				linear.add(-1,"k_"+startid+"_"+i);
				problem.add(new Constraint("", linear, ">=",0));
				
				linear = new Linear();
				linear.add(1,"f_"+i);
				linear.add(1,"k_"+startid+"_"+i);
				problem.add(new Constraint("", linear, "<=",1));
				
				linear = new Linear();
				linear.add(1,"f_"+i);
				linear.add(1,"k_"+startid+"_"+i);
				linear.add(-1,"besr_"+startid);
				problem.add(new Constraint("", linear, ">=",0));
				
				//head_destination = head_i_j - headloss
				linear = new Linear();
				linear.add(1,"headloss_"+i);
				linear.add(1,"head_"+destinationid);
				linear.add(-1,"head_"+startid+"_"+i);
				problem.add(new Constraint("", linear, "=",0));
			}
			else{
				// head_destination = head_source - headloss
				linear = new Linear();
				linear.add(1,"headloss_"+i);
				linear.add(1,"head_"+destinationid);
				linear.add(-1,"head_"+startid);
				problem.add(new Constraint("", linear, "=",0));
			}
		}
	}
	
	//set headloss constraints in pipes
	//includes ESRs and thus primary and secondary network considerations
	//also allows nodes with zero demands to have ESRs
	private void setHeadLossConstraints_esr_gen() throws Exception{
		Linear linear = new Linear();

		for(Node node : nodes.values()){
			int i = node.getNodeID();
			if(source==node){
				linear = new Linear();
				linear.add(1,"head_"+i);
				problem.add(new Constraint("", linear, "=", source.getHead()));
			}
			else{
				linear = new Linear();
				linear.add(1,"head_"+i);
				linear.add(-1,"esr_"+i);
				problem.add(new Constraint("", linear, ">=", node.getElevation() + node.getResidualPressure()));
			
			}
						
			// min_esr_height <= esr_i <= max_esr_height
			linear = new Linear();
			linear.add(minEsrHeight,"besr_"+i);
			linear.add(-1,"esr_"+i);
			problem.add(new Constraint("", linear, "<=",0));
							
			linear = new Linear();
			linear.add(maxEsrHeight,"besr_"+i);
			linear.add(-1,"esr_"+i);
			problem.add(new Constraint("", linear, ">=",0));
			
			linear = new Linear();
			linear.add(1,"besr_"+i);
			linear.add(-1,"s_"+i+"_"+i);
			problem.add(new Constraint("", linear, "=",0));
			
		}		
		
		for(Pipe pipe : pipes.values()){
			linear = new Linear();
			int i = pipe.getPipeID();
			int j = 0;
			if(pipe.isAllowParallel()){
				for(PipeCost entry : pipeCost){
					// primary flow in case of parallel pipe
					double flow = pipe.getFlow() / (1 + (entry.getRoughness()/pipe.getRoughness())*Math.pow(entry.getDiameter()/pipe.getDiameter(), 4.87/1.852));
					double headloss = Util.HWheadLoss(pipe.getLength(), flow, pipe.getRoughness(), pipe.getDiameter());
					double headlossSecondary = Util.HWheadLoss(pipe.getLength(), flow * secondaryFlowFactor, pipe.getRoughness(), pipe.getDiameter());
					linear.add(headlossSecondary, "p_"+i+"_"+j);
					linear.add(headloss-headlossSecondary, "yp_"+i+"_"+j);
					j++;
				}	
				// primary flow in case of no parallel pipe
				double flow = pipe.getFlow(); 
				double headloss = Util.HWheadLoss(pipe.getLength(), flow, pipe.getRoughness(), pipe.getDiameter());
				double headlossSecondary = Util.HWheadLoss(pipe.getLength(), flow * secondaryFlowFactor, pipe.getRoughness(), pipe.getDiameter());
				linear.add(headlossSecondary, "p_"+i);
				linear.add(headloss-headlossSecondary, "yp_"+i);					
			}
			else{
				if(pipe.getDiameter()!=0){
					j=0;
					for(PipeCost entry : pipeCost){
						double headloss = Util.HWheadLoss(pipe.getFlow(), pipe.getRoughness(), entry.getDiameter());
						double headlossSecondary = Util.HWheadLoss(pipe.getFlow() * secondaryFlowFactor, pipe.getRoughness(), entry.getDiameter());
						linear.add(headlossSecondary, "l_"+i+"_"+j);
						linear.add(headloss-headlossSecondary, "y_"+i+"_"+j);
						j++;
					}
				}
				else{	
					j=0;
					for(PipeCost entry : pipeCost){
						double headloss = Util.HWheadLoss(pipe.getFlow(), entry.getRoughness(), entry.getDiameter());
						double headlossSecondary = Util.HWheadLoss(pipe.getFlow() * secondaryFlowFactor, entry.getRoughness(), entry.getDiameter());
						linear.add(headlossSecondary, "l_"+i+"_"+j);
						linear.add(headloss-headlossSecondary, "y_"+i+"_"+j);
						j++;
					}	
				}
			}
			Node startNode = pipe.getStartNode();
			int startid = startNode.getNodeID();
			int destinationid = pipe.getEndNode().getNodeID();
			
			linear.add(-1,"headloss_"+i);
			problem.add(new Constraint("", linear, "=",0));
			
			
			// head_source = head_source if demand_source = 0
			// head_source = head_elevation + esr_height if s_source_source = 1 and f_pipe = 0
			// else head_source = head_source
			// introduce head_i_j : is the start node, j is the id of the pipe
			// k_i_j = s_i_i && !f_j
			// head_i_j = k_i_j * (elevation_i + esr_i - head_i) + head_i
			// ehead_i = elevation_i + esr_i - head_i
			// yhead_i_j = k_i_j * ehead_i
					
			// head_i_j = yhead_i_j + head_i
			linear = new Linear();
			linear.add(1,"yhead_"+startid+"_"+i);
			linear.add(1,"head_"+startid);
			linear.add(-1,"head_"+startid+"_"+i);
			problem.add(new Constraint("", linear, "=",0));
			
			// ehead_i = elevation_i + esr_i - head_i
			linear = new Linear();
			linear.add(1,"ehead_"+startid);
			linear.add(1,"head_"+startid);
			linear.add(-1,"esr_"+startid);
			problem.add(new Constraint("", linear, "=",startNode.getElevation()));
			
			//to capture yhead_startid_i = k_startid_i * ehead_startid
			// ehead range assumed to be -10000,10000
			int minEHead = -10000;
			int maxEHead = 10000;
					
			linear = new Linear();
			linear.add(1,"yhead_"+startid+"_"+i);
			problem.add(new Constraint("", linear, "<=",maxEHead));
			
			linear = new Linear();
			linear.add(1,"yhead_"+startid+"_"+i);
			problem.add(new Constraint("", linear, ">=",minEHead));
			
			linear = new Linear();
			linear.add(-1,"yhead_"+startid+"_"+i);
			linear.add(maxEHead,"k_"+startid+"_"+i);
			problem.add(new Constraint("", linear, ">=",0));
			
			linear = new Linear();
			linear.add(-1,"yhead_"+startid+"_"+i);
			linear.add(minEHead,"k_"+startid+"_"+i);
			problem.add(new Constraint("", linear, "<=",0));
			
			linear = new Linear();
			linear.add(1,"ehead_"+startid);
			linear.add(maxEHead,"k_"+startid+"_"+i);
			linear.add(-1,"yhead_"+startid+"_"+i);
			problem.add(new Constraint("", linear, "<=",maxEHead));
			
			linear = new Linear();
			linear.add(1,"ehead_"+startid);
			linear.add(minEHead,"k_"+startid+"_"+i);
			linear.add(-1,"yhead_"+startid+"_"+i);
			problem.add(new Constraint("", linear, ">=",minEHead));
			
			linear = new Linear();
			linear.add(-1,"ehead_"+startid);
			linear.add(maxEHead,"k_"+startid+"_"+i);
			linear.add(1,"yhead_"+startid+"_"+i);
			problem.add(new Constraint("", linear, "<=",maxEHead));
			
			
			//to capture k_startid_i = s_startid_startid ^ !f_i
			
			linear = new Linear();
			linear.add(1,"s_"+startid+"_"+startid);
			linear.add(-1,"k_"+startid+"_"+i);
			problem.add(new Constraint("", linear, ">=",0));
			
			linear = new Linear();
			linear.add(1,"f_"+i);
			linear.add(1,"k_"+startid+"_"+i);
			problem.add(new Constraint("", linear, "<=",1));
			
			linear = new Linear();
			linear.add(1,"f_"+i);
			linear.add(1,"k_"+startid+"_"+i);
			linear.add(-1,"s_"+startid+"_"+startid);
			problem.add(new Constraint("", linear, ">=",0));
			
			//head_destination = head_i_j - headloss
			linear = new Linear();
			linear.add(1,"headloss_"+i);
			linear.add(1,"head_"+destinationid);
			linear.add(-1,"head_"+startid+"_"+i);
			problem.add(new Constraint("", linear, "=",0));
		}
	}
	
	//add constraints related to the ESR structure in the network
	//ensures ESR configuration is valid
	//includes constraints related to selection of cost row table for each ESR
	private void setEsrConstraints() throws Exception{		
		Linear linear,linear2;
		for(Node node : nodes.values()){
			int i = node.getNodeID();
			if(node.getDemand()!=0){
				
				for(Node downNode : node.getDownstreamNodes()){
					int j = downNode.getNodeID();
					
					//following constraint represents s_i_i = 0 => s_j_j = 0 where j are downstream of i 
					linear = new Linear();
					linear.add(1, "s_"+i+"_"+i); 
					linear.add(-1, "s_"+j+"_"+j);
					problem.add(new Constraint("", linear, ">=",0));
					
					//following constraint represents s_i_i = 0 => s_i_j = 0 where j are downstream of i 
					linear = new Linear();
					linear.add(1, "s_"+i+"_"+i); 
					linear.add(-1, "s_"+i+"_"+j);
					problem.add(new Constraint("", linear, ">=",0));
				}
				
				// following constraint represents sum over parent i's s_j_i = 1 
				// i.e. exactly one esr serves location i 
				linear = new Linear();
				linear.add(1, "s_"+i+"_"+i);
				for(Node upstreamNode : node.getUpstreamNodes()){
					int j = upstreamNode.getNodeID();
					linear.add(1, "s_"+j+"_"+i);
				}
				problem.add(new Constraint("", linear, "=",1));
		
				// d_i = sum all demand served by i / esrcapacity factor
				linear = new Linear();
				linear.add(node.getRequiredCapacity(esrCapacityFactor), "s_"+i+"_"+i);
				for(Node downstreamNode : node.getDownstreamNodes()){
					int j = downstreamNode.getNodeID();
					linear.add(downstreamNode.getRequiredCapacity(esrCapacityFactor), "s_"+i+"_"+j);
				}
				linear.add(-1,"d_"+i);
				problem.add(new Constraint("", linear, "=",0));
				
				int j=0;
				linear2 = new Linear();
				for(EsrCost esr : esrCost){
					
					//to determine e_i_j
					linear = new Linear();
					linear.add(esr.getMinCapacity(),"e_"+i+"_"+j);
					linear.add(-1,"d_"+i);
					problem.add(new Constraint("", linear, "<=",0));
					
					linear = new Linear();
					linear.add(totalDemand,"e_"+i+"_"+j); //totaldemand just stands for a large amount of demand M
					linear.add(1,"d_"+i);
					problem.add(new Constraint("", linear, "<=",esr.getMaxCapacity()+totalDemand));
					
					linear2.add(1,"e_"+i+"_"+j);
					
					//to capture z_i_j = d_i * e_i_j
					
					linear = new Linear();
					linear.add(totalDemand,"e_"+i+"_"+j);
					linear.add(-1,"z_"+i+"_"+j);
					problem.add(new Constraint("", linear, ">=",0));
					
					linear = new Linear();
					linear.add(1,"d_"+i);
					linear.add(-1,"z_"+i+"_"+j);
					problem.add(new Constraint("", linear, ">=",0));
					
					linear = new Linear();
					linear.add(totalDemand,"e_"+i+"_"+j);
					linear.add(-1,"z_"+i+"_"+j);
					linear.add(1,"d_"+i);
					problem.add(new Constraint("", linear, "<=",totalDemand));
					
					j++;
				}
				
				problem.add(new Constraint("", linear2, "=",1));
			}
		}
	}
	
	//add constraints related to the ESR structure in the network
	//ensures ESR configuration is valid
	//includes constraints related to selection of cost row table for each ESR
	//allows zero demand nodes
	private void setEsrConstraints_gen() throws Exception{		
		Linear linear,linear2;
		for(Node node : nodes.values()){
			int i = node.getNodeID();
			
			for(Pipe pipe : node.getOutgoingPipes()){
				
				int j = pipe.getEndNode().getNodeID();
				
				//following constraint represents s_i_i = 0 => s_j_j = 0 where j are children of i 
				linear = new Linear();
				linear.add(1, "s_"+i+"_"+i); 
				linear.add(-1, "s_"+j+"_"+j);
				problem.add(new Constraint("", linear, ">=",0));
			}
			
			// following constraint represents sum over parent i's s_j_i = 1 
			// i.e. exactly one esr serves location i 
			linear = new Linear();
			linear.add(1, "s_"+i+"_"+i);
			for(Node upstreamNode : node.getUpstreamNodes()){
				int j = upstreamNode.getNodeID();
				linear.add(1, "s_"+j+"_"+i);
			}
			problem.add(new Constraint("", linear, "=",1));
	
			// d_i = sum all demand served by i / esrcapacity factor
			linear = new Linear();
			linear.add(node.getRequiredCapacity(esrCapacityFactor), "s_"+i+"_"+i);
			for(Node downstreamNode : node.getDownstreamNodes()){
				int j = downstreamNode.getNodeID();
				linear.add(downstreamNode.getRequiredCapacity(esrCapacityFactor), "s_"+i+"_"+j);
			}
			linear.add(-1,"d_"+i);
			problem.add(new Constraint("", linear, "=",0));
			
			int j=0;
			linear2 = new Linear();
			for(EsrCost esr : esrCost){
				
				//to determine e_i_j
				linear = new Linear();
				linear.add(esr.getMinCapacity(),"e_"+i+"_"+j);
				linear.add(-1,"d_"+i);
				problem.add(new Constraint("", linear, "<=",0));
				
				linear = new Linear();
				linear.add(totalDemand,"e_"+i+"_"+j); //totaldemand just stands for a large amount of demand M
				linear.add(1,"d_"+i);
				problem.add(new Constraint("", linear, "<=",esr.getMaxCapacity()+totalDemand));
				
				linear2.add(1,"e_"+i+"_"+j);
				
				//to capture z_i_j = d_i * e_i_j
				
				linear = new Linear();
				linear.add(totalDemand,"e_"+i+"_"+j);
				linear.add(-1,"z_"+i+"_"+j);
				problem.add(new Constraint("", linear, ">=",0));
				
				linear = new Linear();
				linear.add(1,"d_"+i);
				linear.add(-1,"z_"+i+"_"+j);
				problem.add(new Constraint("", linear, ">=",0));
				
				linear = new Linear();
				linear.add(totalDemand,"e_"+i+"_"+j);
				linear.add(-1,"z_"+i+"_"+j);
				linear.add(1,"d_"+i);
				problem.add(new Constraint("", linear, "<=",totalDemand));
				
				j++;
			}
			
			problem.add(new Constraint("", linear2, "=",1));
		}
		
		for(Pipe pipe : pipes.values()){
			
			Node parent = pipe.getStartNode();
			Node node = pipe.getEndNode();
			
			int i = node.getNodeID();
			int p = parent.getNodeID();
			
			// represents s_i_i = 0 => (s_k_p = s_k_i) where p is the parent of i
			linear = new Linear();
			linear2 = new Linear();
			int M = 0;
			int counter = 1;
			
			Set<Node> candidates = new HashSet<Node>();
			candidates.addAll(parent.getUpstreamNodes());
			candidates.add(parent);
			
			for(Node n : candidates){
				int k = n.getNodeID();
									
				linear.add(counter, "s_"+k+"_"+i);
				linear.add(-1 * counter, "s_"+k+"_"+p);
				
				linear2.add(counter, "s_"+k+"_"+i);
				linear2.add(-1 * counter, "s_"+k+"_"+p);
				
				M = M + counter;
				counter = counter * 2;
			}
			
			linear.add(-1*M,"s_"+i+"_"+i);
			linear2.add(M,"s_"+i+"_"+i);
			
			problem.add(new Constraint("", linear, "<=",0));			
			problem.add(new Constraint("", linear2, ">=",0));
		}
		
	}
	
	//add constraints related to the ESR structure in the network
	//ensures ESR configuration is valid
	//includes constraints related to selection of cost row table for each ESR
	//allows removal of certain nodes from ESR consideration
	private void setEsrConstraints_gen2() throws Exception{		
		Linear linear,linear2;
		for(Node node : nodes.values()){
			int i = node.getNodeID();
			
			for(Pipe pipe : node.getOutgoingPipes()){
				
				int j = pipe.getEndNode().getNodeID();
				
				//following constraint represents s_i_i = 0 => s_j_j = 0 where j are children of i 
				linear = new Linear();
				linear.add(1, "s_"+i+"_"+i); 
				linear.add(-1, "s_"+j+"_"+j);
				problem.add(new Constraint("", linear, ">=",0));
			}
			
			// following constraint represents sum over parent i's s_j_i = 1 
			// i.e. exactly one esr serves location i 
			linear = new Linear();
			linear.add(1, "s_"+i+"_"+i);
			for(Node upstreamNode : node.getUpstreamNodes()){
				int j = upstreamNode.getNodeID();
				linear.add(1, "s_"+j+"_"+i);
			}
			problem.add(new Constraint("", linear, "=",1));
	
			// d_i = sum all demand served by i / esrcapacity factor
			linear = new Linear();
			linear.add(node.getRequiredCapacity(esrCapacityFactor), "s_"+i+"_"+i);
			for(Node downstreamNode : node.getDownstreamNodes()){
				int j = downstreamNode.getNodeID();
				linear.add(downstreamNode.getRequiredCapacity(esrCapacityFactor), "s_"+i+"_"+j);
			}
			linear.add(-1,"d_"+i);
			problem.add(new Constraint("", linear, "=",0));
			
			
			if(node.getAllowESR()){
				int j=0;
				linear2 = new Linear();
				for(EsrCost esr : esrCost){
					
					//to determine e_i_j
					linear = new Linear();
					linear.add(esr.getMinCapacity(),"e_"+i+"_"+j);
					linear.add(-1,"d_"+i);
					problem.add(new Constraint("", linear, "<=",0));
					
					linear = new Linear();
					linear.add(totalDemand,"e_"+i+"_"+j); //totaldemand just stands for a large amount of demand M
					linear.add(1,"d_"+i);
					problem.add(new Constraint("", linear, "<=",esr.getMaxCapacity()+totalDemand));
					
					linear2.add(1,"e_"+i+"_"+j);
					
					//to capture z_i_j = d_i * e_i_j
					
					linear = new Linear();
					linear.add(totalDemand,"e_"+i+"_"+j);
					linear.add(-1,"z_"+i+"_"+j);
					problem.add(new Constraint("", linear, ">=",0));
					
					linear = new Linear();
					linear.add(1,"d_"+i);
					linear.add(-1,"z_"+i+"_"+j);
					problem.add(new Constraint("", linear, ">=",0));
					
					linear = new Linear();
					linear.add(totalDemand,"e_"+i+"_"+j);
					linear.add(-1,"z_"+i+"_"+j);
					linear.add(1,"d_"+i);
					problem.add(new Constraint("", linear, "<=",totalDemand));
					
					j++;
				}
				
				problem.add(new Constraint("", linear2, "=",1));
			}
		}
		
		for(Pipe pipe : pipes.values()){
			
			Node parent = pipe.getStartNode();
			Node node = pipe.getEndNode();
			
			int i = node.getNodeID();
			int p = parent.getNodeID();
			
			// represents s_i_i = 0 => (s_k_p = s_k_i) where p is the parent of i
			linear = new Linear();
			linear2 = new Linear();
			int M = 0;
			int counter = 1;
			
			Set<Node> candidates = new HashSet<Node>();
			
			for(Node n : node.getUpstreamNodes()){
				if(n.getAllowESR())
					candidates.add(n);
			}
			
			for(Node n : candidates){
				int k = n.getNodeID();
									
				linear.add(counter, "s_"+k+"_"+i);
				linear.add(-1 * counter, "s_"+k+"_"+p);
				
				linear2.add(counter, "s_"+k+"_"+i);
				linear2.add(-1 * counter, "s_"+k+"_"+p);
				
				M = M + counter;
				counter = counter * 2;
			}
			
			linear.add(-1*M,"s_"+i+"_"+i);
			linear2.add(M,"s_"+i+"_"+i);
			
			problem.add(new Constraint("", linear, "<=",0));			
			problem.add(new Constraint("", linear2, ">=",0));
		}
		
	}
		
	//add constraints related to the ESR structure in the network
	//ensures ESR configuration is valid
	//includes constraints related to selection of cost row table for each ESR
	//alternative z_i_j constraints
	private void setEsrConstraints_gen4() throws Exception{		
		Linear linear,linear2;
		for(Node node : nodes.values()){
			int i = node.getNodeID();
			
			for(Pipe pipe : node.getOutgoingPipes()){
				
				int j = pipe.getEndNode().getNodeID();
				
				//following constraint represents s_i_i = 0 => s_j_j = 0 where j are children of i 
				linear = new Linear();
				linear.add(1, "s_"+i+"_"+i); 
				linear.add(-1, "s_"+j+"_"+j);
				problem.add(new Constraint("", linear, ">=",0));
			}
			
			// following constraint represents sum over parent i's s_j_i = 1 
			// i.e. exactly one esr serves location i 
			linear = new Linear();
			linear.add(1, "s_"+i+"_"+i);
			for(Node upstreamNode : node.getUpstreamNodes()){
				int j = upstreamNode.getNodeID();
				linear.add(1, "s_"+j+"_"+i);
			}
			problem.add(new Constraint("", linear, "=",1));
	
			// d_i = sum all demand served by i / esrcapacity factor
			linear = new Linear();
			linear.add(node.getRequiredCapacity(esrCapacityFactor), "s_"+i+"_"+i);
			for(Node downstreamNode : node.getDownstreamNodes()){
				int j = downstreamNode.getNodeID();
				linear.add(downstreamNode.getRequiredCapacity(esrCapacityFactor), "s_"+i+"_"+j);
			}
			linear.add(-1,"d_"+i);
			problem.add(new Constraint("", linear, "=",0));
			
			
			if(node.getAllowESR()){
				int j=0;
				linear2 = new Linear();
				Linear linear3 = new Linear();
				for(EsrCost esr : esrCost){
					
					//to determine e_i_j and z_i_j
					// e_i_j*min_j <= z_i_j <= e_i_j*max_j
					// sum e_i_j = 1
					// sum z_i_j = d_i
					
					linear = new Linear();
					linear.add(esr.getMinCapacity(),"e_"+i+"_"+j);
					linear.add(-1,"z_"+i+"_"+j);
					problem.add(new Constraint("", linear, "<=",0));
					
					linear = new Linear();
					linear.add(esr.getMaxCapacity(),"e_"+i+"_"+j);
					linear.add(-1,"z_"+i+"_"+j);
					problem.add(new Constraint("", linear, ">=",0));
										
					linear2.add(1,"e_"+i+"_"+j);
					linear3.add(1,"z_"+i+"_"+j);
					
					j++;
				}
				
				problem.add(new Constraint("", linear2, "=",1));
				
				linear3.add(-1,"d_"+i);
				problem.add(new Constraint("", linear3, "=",0));
			}
		}
		
		for(Pipe pipe : pipes.values()){
			
			Node parent = pipe.getStartNode();
			Node node = pipe.getEndNode();
			
			int i = node.getNodeID();
			int p = parent.getNodeID();
			
			// represents s_i_i = 0 => (s_k_p = s_k_i) where p is the parent of i
			linear = new Linear();
			linear2 = new Linear();
			int M = 0;
			int counter = 1;
			
			Set<Node> candidates = new HashSet<Node>();
			
			for(Node n : node.getUpstreamNodes()){
				if(n.getAllowESR())
					candidates.add(n);
			}
			
			for(Node n : candidates){
				int k = n.getNodeID();
									
				linear.add(counter, "s_"+k+"_"+i);
				linear.add(-1 * counter, "s_"+k+"_"+p);
				
				linear2.add(counter, "s_"+k+"_"+i);
				linear2.add(-1 * counter, "s_"+k+"_"+p);
				
				M = M + counter;
				counter = counter * 2;
			}
			
			linear.add(-1*M,"s_"+i+"_"+i);
			linear2.add(M,"s_"+i+"_"+i);
			
			problem.add(new Constraint("", linear, "<=",0));			
			problem.add(new Constraint("", linear2, ">=",0));
		}
		
	}
	
	//add constraints related to the ESR structure in the network
	//ensures ESR configuration is valid
	//includes constraints related to selection of cost row table for each ESR
	//replace how to implement constraint s_i_i = 0 => s_k_i = s_k_p
	private void setEsrConstraints_gen5() throws Exception{		
		Linear linear,linear2;
		for(Node node : nodes.values()){
			int i = node.getNodeID();
			
			for(Pipe pipe : node.getOutgoingPipes()){
				
				int j = pipe.getEndNode().getNodeID();
				
				//following constraint represents s_i_i = 0 => s_j_j = 0 where j are children of i 
				linear = new Linear();
				linear.add(1, "s_"+i+"_"+i); 
				linear.add(-1, "s_"+j+"_"+j);
				problem.add(new Constraint("", linear, ">=",0));
			}
			
			// following constraint represents sum over parent i's s_j_i = 1 
			// i.e. exactly one esr serves location i 
			linear = new Linear();
			linear.add(1, "s_"+i+"_"+i);
			for(Node upstreamNode : node.getUpstreamNodes()){
				int j = upstreamNode.getNodeID();
				linear.add(1, "s_"+j+"_"+i);
			}
			problem.add(new Constraint("", linear, "=",1));
	
			// d_i = sum all demand served by i / esrcapacity factor
			linear = new Linear();
			linear.add(node.getRequiredCapacity(esrCapacityFactor), "s_"+i+"_"+i);
			for(Node downstreamNode : node.getDownstreamNodes()){
				int j = downstreamNode.getNodeID();
				linear.add(downstreamNode.getRequiredCapacity(esrCapacityFactor), "s_"+i+"_"+j);
			}
			linear.add(-1,"d_"+i);
			problem.add(new Constraint("", linear, "=",0));
			
			
			if(node.getAllowESR()){
				int j=0;
				linear2 = new Linear();
				Linear linear3 = new Linear();
				for(EsrCost esr : esrCost){
					
					//to determine e_i_j and z_i_j
					// e_i_j*min_j <= z_i_j <= e_i_j*max_j
					// sum e_i_j = 1
					// sum z_i_j = d_i
					
					linear = new Linear();
					linear.add(esr.getMinCapacity(),"e_"+i+"_"+j);
					linear.add(-1,"z_"+i+"_"+j);
					problem.add(new Constraint("", linear, "<=",0));
					
					linear = new Linear();
					linear.add(esr.getMaxCapacity(),"e_"+i+"_"+j);
					linear.add(-1,"z_"+i+"_"+j);
					problem.add(new Constraint("", linear, ">=",0));
										
					linear2.add(1,"e_"+i+"_"+j);
					linear3.add(1,"z_"+i+"_"+j);
					
					j++;
				}
				
				problem.add(new Constraint("", linear2, "=",1));
				
				linear3.add(-1,"d_"+i);
				problem.add(new Constraint("", linear3, "=",0));
			}
		}
		
		for(Pipe pipe : pipes.values()){
			
			Node parent = pipe.getStartNode();
			Node node = pipe.getEndNode();
			
			int i = node.getNodeID();
			int p = parent.getNodeID();
			
			// represents s_i_i = 0 => (s_k_p = s_k_i) where p is the parent of i			
			Set<Node> candidates = new HashSet<Node>();
			
			for(Node n : node.getUpstreamNodes()){
				if(n.getAllowESR())
					candidates.add(n);
			}
			
			for(Node n : candidates){
				int k = n.getNodeID();
				
				linear = new Linear();
				linear.add(1, "s_"+k+"_"+i);
				linear.add(-1, "s_"+k+"_"+p);
				linear.add(-1, "s_"+i+"_"+i);
				problem.add(new Constraint("", linear, "<=",0));
				
				linear = new Linear();
				linear.add(-1, "s_"+k+"_"+i);
				linear.add(1, "s_"+k+"_"+p);
				linear.add(-1, "s_"+i+"_"+i);
				problem.add(new Constraint("", linear, "<=",0));
			}
		}
	}
	
	//add constraints related to user defined list of nodes that must have ESR / cannot have ESR
	private void setEsrOptionConstraints() throws Exception{
		if(!esrGeneralProperties.allow_dummy){
			for(Node n : nodes.values()){
				if(n.getDemand()==0)
					n.setAllowESR(false);
			}
		}
		
		if(esrGeneralProperties.must_esr!=null){
			for(int nodeid : esrGeneralProperties.must_esr){
				Node node = nodes.get(nodeid);
				if(node==null){
	            	throw new Exception("Invalid node:" + nodeid + " provided for must have ESR list in ESR options");
	            }
				Linear linear = new Linear();
				linear.add(1, "s_"+nodeid+"_"+nodeid); 
				problem.add(new Constraint("", linear, "=",1));
				
				linear = new Linear();
				linear.add(1, "d_"+nodeid); 
				problem.add(new Constraint("", linear, ">=",0.00001));
			}
		}
		
		if(esrGeneralProperties.must_not_esr!=null){
			for(int nodeid : esrGeneralProperties.must_not_esr){
				Node node = nodes.get(nodeid);
				if(node==null){
	            	throw new Exception("Invalid node:" + nodeid + " provided for must not have ESR list in ESR options");
	            }
				node.setAllowESR(false);
			}
		}
	}

	//add constraints related to user defined list of nodes that must have ESR / cannot have ESR
	private void setEsrOptionConstraints_gen2() throws Exception{
		if(!esrGeneralProperties.allow_dummy){
			for(Node n : nodes.values()){
				if(n.getDemand()==0)
					n.setAllowESR(false);
			}
		}
		
		if(esrGeneralProperties.must_esr!=null){
			for(int nodeid : esrGeneralProperties.must_esr){
				Node node = nodes.get(nodeid);
				if(node==null){
	            	throw new Exception("Invalid node:" + nodeid + " provided for must have ESR list in ESR options");
	            }
				Linear linear = new Linear();
				
				linear = new Linear();
				linear.add(1, "d_"+nodeid); 
				problem.add(new Constraint("", linear, ">=",0.00001));
			}
		}
		
		if(esrGeneralProperties.must_not_esr!=null){
			for(int nodeid : esrGeneralProperties.must_not_esr){
				Node node = nodes.get(nodeid);
				if(node==null){
	            	throw new Exception("Invalid node:" + nodeid + " provided for must not have ESR list in ESR options");
	            }
				node.setAllowESR(false);
			}
		}
	}
	
	//set valve settings in pipes
	private void setValveOptionConstraints() throws Exception{
		if(valves!=null){
			for(ValveStruct valve:valves){
				Pipe pipe = pipes.get(valve.pipeid);
				if(pipe==null){
	            	throw new Exception("Invalid pipe:" + valve.pipeid + " provided in pressure reducing valve list");
	            }
				pipe.setValveSetting(valve.valvesetting);
			}
		}
	}
	
	//add pump manual option related constraints
	private void setPumpOptionConstraints() throws Exception{
		if(pumpGeneralProperties.pump_enabled){
			if(pumpGeneralProperties.must_not_pump!=null){
				for(int pipeid : pumpGeneralProperties.must_not_pump){
					Pipe pipe = pipes.get(pipeid);
					if(pipe==null){
		            	throw new Exception("Invalid pipe:" + pipeid + " provided for must not have pump list in pump options");
		            }
					pipe.setAllowPump(false);
				}
			}
			
			for(PumpManualStruct s: pumpManualArray){
				Pipe pipe = pipes.get(s.pipeid);
				if(pipe==null){
	            	throw new Exception("Invalid pipe:" + s.pipeid + " provided for manual pump list in pump options");
	            }
				Linear linear = new Linear();
				linear.add(1, "pumppower_"+s.pipeid);
				problem.add(new Constraint("", linear, "=", s.pumppower));
			}
		}
	}
	
	//add constraints related to the ESR structure in the network
	//ensures ESR configuration is valid
	//includes constraints related to selection of cost row table for each ESR
	//replace how to implement constraint s_i_i = 0 => s_k_i = s_k_p
	private void setEsrConstraints_gen6() throws Exception{		
		Linear linear,linear2;
		for(Node node : nodes.values()){
			int i = node.getNodeID();
			
			for(Pipe pipe : node.getOutgoingPipes()){
				
				int j = pipe.getEndNode().getNodeID();
				
				//following constraint represents s_i_i = 0 => s_j_j = 0 where j are children of i 
				linear = new Linear();
				linear.add(1, "s_"+i+"_"+i); 
				linear.add(-1, "s_"+j+"_"+j);
				problem.add(new Constraint("", linear, ">=",0));
				
				// along a certain pipe "pipeid", sum(s_i_j) >=1 => f_pipeid = 0
				int n = 1;
				linear = new Linear();
				linear.add(1,"s_"+i+"_"+j);
				for(Node downNode : pipe.getEndNode().getDownstreamNodes()){
					int k = downNode.getNodeID();
					linear.add(1,"s_"+i+"_"+k);
					n++;
				}
				linear.add(n,"f_"+pipe.getPipeID());
				problem.add(new Constraint("", linear, "<=",n));
				
			}
			
			// following constraint represents sum over parent i's s_j_i = 1 
			// i.e. exactly one esr serves location i 
			linear = new Linear();
			linear.add(1, "s_"+i+"_"+i);
			for(Node upstreamNode : node.getUpstreamNodes()){
				int j = upstreamNode.getNodeID();
				linear.add(1, "s_"+j+"_"+i);
			}
			problem.add(new Constraint("", linear, "=",1));
	
			// d_i = sum all demand served by i / esrcapacity factor
			double totaldemand = 0;
			totaldemand = totaldemand + node.getRequiredCapacity(esrCapacityFactor);
			linear = new Linear();
			linear.add(node.getRequiredCapacity(esrCapacityFactor), "s_"+i+"_"+i);
			for(Node downstreamNode : node.getDownstreamNodes()){
				int j = downstreamNode.getNodeID();
				totaldemand = totaldemand + downstreamNode.getRequiredCapacity(esrCapacityFactor);
				linear.add(downstreamNode.getRequiredCapacity(esrCapacityFactor), "s_"+i+"_"+j);
			}
			linear.add(-1,"d_"+i);
			//problem.add(new Constraint("", linear, "=",0));
			problem.add(new Constraint("", linear, "<=",0.0001));
			problem.add(new Constraint("", linear, ">=",-0.0001));
			
			// s_i_i = 0 => d_i=0 i.e. all s_i_j=0
			linear = new Linear();
			linear.add(-1,"d_"+i);
			linear.add(totaldemand,"s_"+i+"_"+i);
			//problem.add(new Constraint("", linear, ">=",0));
			problem.add(new Constraint("", linear, ">=",-0.0001));
			
			if(node.getAllowESR()){
				int j=0;
				linear2 = new Linear();
				Linear linear3 = new Linear();
				for(EsrCost esr : esrCost){
					
					//to determine e_i_j and z_i_j
					// e_i_j*min_j <= z_i_j <= e_i_j*max_j
					// sum e_i_j = 1
					// sum z_i_j = d_i
					
					linear = new Linear();
					linear.add(esr.getMinCapacity(),"e_"+i+"_"+j);
					linear.add(-1,"z_"+i+"_"+j);
					problem.add(new Constraint("", linear, "<=",0));
					
					linear = new Linear();
					linear.add(esr.getMaxCapacity(),"e_"+i+"_"+j);
					linear.add(-1,"z_"+i+"_"+j);
					problem.add(new Constraint("", linear, ">=",0));
										
					linear2.add(1,"e_"+i+"_"+j);
					linear3.add(1,"z_"+i+"_"+j);
					
					j++;
				}
				
				problem.add(new Constraint("", linear2, "=",1));
				
				linear3.add(-1,"d_"+i);
				problem.add(new Constraint("", linear3, "=",0));
			}
		}		
	}

	//add constraints related to the ESR structure in the network
	//ensures ESR configuration is valid
	//includes constraints related to selection of cost row table for each ESR
	//replace how to implement constraint s_i_i = 0 => s_k_i = s_k_p
	private void setEsrConstraints_gen7() throws Exception{		
		Linear linear,linear2;
		for(Node node : nodes.values()){
			int i = node.getNodeID();
			
			for(Pipe pipe : node.getOutgoingPipes()){
				
				int j = pipe.getEndNode().getNodeID();
				
				//following constraint represents s_i_i = 0 => s_j_j = 0 where j are children of i 
				linear = new Linear();
				linear.add(1, "s_"+i+"_"+i); 
				linear.add(-1, "s_"+j+"_"+j);
				problem.add(new Constraint("", linear, ">=",0));
				
				
				//following constraint represents s_i_i = 0 => s_i_j = 0 where j are children of i 
				//linear = new Linear();
				//linear.add(1, "s_"+i+"_"+i); 
				//linear.add(-1, "s_"+i+"_"+j);
				//problem.add(new Constraint("", linear, ">=",0));
				
				
				// following constraint represents s_i_j = s_i_k where j is child of i and k are downstream of j
				for(Node downNode : pipe.getEndNode().getDownstreamNodes()){
					linear = new Linear();
					linear.add(1,"s_"+i+"_"+j);
					
					int k = downNode.getNodeID();
					linear.add(-1,"s_"+i+"_"+k);
					problem.add(new Constraint("", linear, "=",0));
				}
				
				
			}
			
			// following constraint represents sum over parent i's s_j_i = 1 
			// i.e. exactly one esr serves location i 
			linear = new Linear();
			linear.add(1, "s_"+i+"_"+i);
			for(Node upstreamNode : node.getUpstreamNodes()){
				int j = upstreamNode.getNodeID();
				linear.add(1, "s_"+j+"_"+i);
			}
			problem.add(new Constraint("", linear, "=",1));
	
			// d_i = sum all demand served by i / esrcapacity factor
			double totaldemand = 0;
			totaldemand = totaldemand + node.getRequiredCapacity(esrCapacityFactor);
			linear = new Linear();
			linear.add(node.getRequiredCapacity(esrCapacityFactor), "s_"+i+"_"+i);
			for(Node downstreamNode : node.getDownstreamNodes()){
				int j = downstreamNode.getNodeID();
				totaldemand = totaldemand + downstreamNode.getRequiredCapacity(esrCapacityFactor);
				linear.add(downstreamNode.getRequiredCapacity(esrCapacityFactor), "s_"+i+"_"+j);
			}
			linear.add(-1,"d_"+i);
			problem.add(new Constraint("", linear, "=",0));
			//problem.add(new Constraint("", linear, "<=",0.0001));
			//problem.add(new Constraint("", linear, ">=",-0.0001));
			
			
			if(node.getAllowESR()){
				int j=0;
				linear2 = new Linear();
				Linear linear3 = new Linear();
				for(EsrCost esr : esrCost){
					
					//to determine e_i_j and z_i_j
					// e_i_j*min_j <= z_i_j <= e_i_j*max_j
					// sum e_i_j = 1
					// sum z_i_j = d_i
					
					linear = new Linear();
					linear.add(esr.getMinCapacity(),"e_"+i+"_"+j);
					linear.add(-1,"z_"+i+"_"+j);
					problem.add(new Constraint("", linear, "<=",0));
					
					linear = new Linear();
					linear.add(esr.getMaxCapacity(),"e_"+i+"_"+j);
					linear.add(-1,"z_"+i+"_"+j);
					problem.add(new Constraint("", linear, ">=",0));
										
					linear2.add(1,"e_"+i+"_"+j);
					linear3.add(1,"z_"+i+"_"+j);
					
					j++;
				}
				
				problem.add(new Constraint("", linear2, "=",1));
				
				linear3.add(-1,"d_"+i);
				problem.add(new Constraint("", linear3, "=",0));
			}
		}		
	}
	
	//add constraints related to the ESR structure in the network
	//ensures ESR configuration is valid
	//includes constraints related to selection of cost row table for each ESR
	//replace how to implement constraint s_i_i = 0 => s_k_i = s_k_p
	private void setEsrConstraints_gen8() throws Exception{		
		Linear linear,linear2;
		for(Node node : nodes.values()){
			int i = node.getNodeID();
			
			linear2 = new Linear();
			linear2.add(node.getRequiredCapacity(esrCapacityFactor), "s_"+i+"_"+i);
			
			for(Pipe pipe : node.getOutgoingPipes()){
				
				int j = pipe.getEndNode().getNodeID();
				
				//d_i sum
				double D_j = 0;
				D_j = pipe.getEndNode().getRequiredCapacity(esrCapacityFactor);
				for(Node n: pipe.getEndNode().getDownstreamNodes()){
					D_j += n.getRequiredCapacity(esrCapacityFactor);
				}
				linear2.add(D_j, "s_"+i+"_"+j);
				
				
				//following constraint represents s_i_i = 0 => s_j_j = 0 where j are children of i 
				//linear = new Linear();
				//linear.add(1, "s_"+i+"_"+i); 
				//linear.add(-1, "s_"+j+"_"+j);
				//problem.add(new Constraint("", linear, ">=",0));
				
				
				//following constraint represents s_i_i = 0 => s_i_j = 0 where j are children of i 
				//linear = new Linear();
				//linear.add(1, "s_"+i+"_"+i); 
				//linear.add(-1, "s_"+i+"_"+j);
				//problem.add(new Constraint("", linear, ">=",0));
				
				
				// following constraint represents s_i_j = s_i_k where j is child of i and k are downstream of j
				for(Node downNode : pipe.getEndNode().getDownstreamNodes()){
					linear = new Linear();
					linear.add(1,"s_"+i+"_"+j);
					
					int k = downNode.getNodeID();
					linear.add(-1,"s_"+i+"_"+k);
					problem.add(new Constraint("", linear, "=",0));
				}
			}
			linear2.add(-1,"d_"+i);
			problem.add(new Constraint("", linear2, "=",0));
			//problem.add(new Constraint("", linear2, "<=",0.0001));
			//problem.add(new Constraint("", linear2, ">=",-0.0001));
			
			
			// following constraint represents sum over parent i's s_j_i = 1 
			// i.e. exactly one esr serves location i 
			linear = new Linear();
			linear.add(1, "s_"+i+"_"+i);
			for(Node upstreamNode : node.getUpstreamNodes()){
				int j = upstreamNode.getNodeID();
				linear.add(1, "s_"+j+"_"+i);
			}
			problem.add(new Constraint("", linear, "=",1));
			
			if(node.getAllowESR()){
				int j=0;
				linear2 = new Linear();
				Linear linear3 = new Linear();
				for(EsrCost esr : esrCost){
					
					//to determine e_i_j and z_i_j
					// e_i_j*min_j <= z_i_j <= e_i_j*max_j
					// sum e_i_j = 1
					// sum z_i_j = d_i
					
					linear = new Linear();
					linear.add(esr.getMinCapacity(),"e_"+i+"_"+j);
					linear.add(-1,"z_"+i+"_"+j);
					problem.add(new Constraint("", linear, "<=",0));
					
					linear = new Linear();
					linear.add(esr.getMaxCapacity(),"e_"+i+"_"+j);
					linear.add(-1,"z_"+i+"_"+j);
					problem.add(new Constraint("", linear, ">=",0));
										
					linear2.add(1,"e_"+i+"_"+j);
					linear3.add(1,"z_"+i+"_"+j);
					
					j++;
				}
				
				problem.add(new Constraint("", linear2, "=",1));
				
				linear3.add(-1,"d_"+i);
				problem.add(new Constraint("", linear3, "=",0));
			}
		}	
	}
	
	//add constraints related to the ESR structure in the network
	//ensures ESR configuration is valid
	//includes constraints related to selection of cost row table for each ESR
	//add constraints to prune ESR cost table rows
	private void setEsrConstraints_gen9() throws Exception{		
		Linear linear,linear2;
		
		for(Node node : nodes.values()){
			int i = node.getNodeID();
			int npipes = node.getOutgoingPipes().size();
			
			
			if(node.getAllowESR()){
				double demand = node.getRequiredCapacity(esrCapacityFactor);
				for(Node n : node.getDownstreamNodes()){
					demand += n.getRequiredCapacity(esrCapacityFactor); 
				}
				int l = 0;
				for(EsrCost e: esrCost){
					if(e.getMinCapacity() > demand){
						linear = new Linear();
						linear.add(1, "e_"+i+"_"+l);
						problem.add(new Constraint("", linear, "=",0));
					}
					l++;
				}
			}
			
			if(node.getAllowESR() && npipes <=0){
				
				for(int k=0;k<Math.pow(2, npipes);k++){
					int mask = 1;
					linear = new Linear();
					linear.add(1, "s_"+i+"_"+i);
					double demand = node.getRequiredCapacity(esrCapacityFactor);
					int rhs = npipes;
					for(Pipe pipe : node.getOutgoingPipes()){
						int j = pipe.getEndNode().getNodeID();
						
						if((mask & k) != 0){	// is current pipe 'chosen'
							linear.add(1, "s_"+i+"_"+j);
							demand += pipe.getEndNode().getRequiredCapacity(esrCapacityFactor);
							for(Node n: pipe.getEndNode().getDownstreamNodes()){
								demand += n.getRequiredCapacity(esrCapacityFactor);
							}
						}
						else{
							linear.add(-1, "s_"+i+"_"+j);
							rhs--;
						}
						mask = mask << 1;
					}
					int l = 0;
					for(EsrCost e: esrCost){
						if(e.getMinCapacity() <= demand && e.getMaxCapacity() > demand){
							linear.add(-1, "e_"+i+"_"+l);
							problem.add(new Constraint("", linear, "<=",rhs));
							//System.out.println(new Constraint("", linear, "<=",rhs));
						}
						l++;
					}
				}
			}
			
			
			if(node.getAllowESR() && npipes <=0){
				
				for(int k=0;k<Math.pow(2, npipes);k++){
					int mask = 1;
					linear = new Linear();
					double demand = node.getRequiredCapacity(esrCapacityFactor);
					int rhs = npipes;
					for(Pipe pipe : node.getOutgoingPipes()){
						if((mask & k) != 0){	// is current pipe 'chosen'
							demand += pipe.getEndNode().getRequiredCapacity(esrCapacityFactor);
							for(Node n: pipe.getEndNode().getDownstreamNodes()){
								demand += n.getRequiredCapacity(esrCapacityFactor);
							}
						}
						mask = mask << 1;
					}
					
					mask = 1;
					linear.add(demand, "s_"+i+"_"+i);
					for(Pipe pipe : node.getOutgoingPipes()){
						int j = pipe.getEndNode().getNodeID();	
						if((mask & k) != 0){	// is current pipe 'chosen'
							linear.add(demand, "s_"+i+"_"+j);
						}
						else{
							linear.add(-1 * demand, "s_"+i+"_"+j);
							rhs--;
						}
						mask = mask << 1;
					}
										
					int l = 0;
					for(EsrCost e: esrCost){
						if(e.getMinCapacity() <= demand && e.getMaxCapacity() > demand){
							linear.add(-1, "z_"+i+"_"+l);
							problem.add(new Constraint("", linear, "<=",rhs*demand));
							//System.out.println(new Constraint("", linear, "<=",rhs));
							break;
						}
						l++;
					}
				}
			}
			
						
			linear2 = new Linear();
			linear2.add(node.getRequiredCapacity(esrCapacityFactor), "s_"+i+"_"+i);
						
			for(Pipe pipe : node.getOutgoingPipes()){
				
				int j = pipe.getEndNode().getNodeID();
								
				//d_i sum
				double D_j = 0;
				D_j = pipe.getEndNode().getRequiredCapacity(esrCapacityFactor);
				for(Node n: pipe.getEndNode().getDownstreamNodes()){
					D_j += n.getRequiredCapacity(esrCapacityFactor);
				}
				linear2.add(D_j, "s_"+i+"_"+j);
				
				//following constraint represents s_i_i = 0 => s_j_j = 0 where j are children of i 
				linear = new Linear();
				linear.add(1, "s_"+i+"_"+i); 
				linear.add(-1, "s_"+j+"_"+j);
				problem.add(new Constraint("", linear, ">=",0));
				
				
				//following constraint represents s_i_i = 0 => s_i_j = 0 where j are children of i 
				linear = new Linear();
				linear.add(1, "s_"+i+"_"+i); 
				linear.add(-1, "s_"+i+"_"+j);
				problem.add(new Constraint("", linear, ">=",0));
				
				
				
				// following constraint represents s_i_j = s_i_k where j is child of i and k are downstream of j
				for(Node downNode : pipe.getEndNode().getDownstreamNodes()){
					linear = new Linear();
					linear.add(1,"s_"+i+"_"+j);
					
					int k = downNode.getNodeID();
					linear.add(-1,"s_"+i+"_"+k);
					problem.add(new Constraint("", linear, "=",0));
				}
			}
			linear2.add(-1,"d_"+i);
			problem.add(new Constraint("", linear2, "=",0));
						
			// following constraint represents sum over parent i's s_j_i = 1 
			// i.e. exactly one esr serves location i 
			linear = new Linear();
			linear.add(1, "s_"+i+"_"+i);
			for(Node upstreamNode : node.getUpstreamNodes()){
				int j = upstreamNode.getNodeID();
				linear.add(1, "s_"+j+"_"+i);
			}
			problem.add(new Constraint("", linear, "=",1));
			
			if(node.getAllowESR()){
				int j=0;
				linear2 = new Linear();
				Linear linear3 = new Linear();
				for(EsrCost esr : esrCost){
					
					//to determine e_i_j and z_i_j
					// e_i_j*min_j <= z_i_j <= e_i_j*max_j
					// sum e_i_j = 1
					// sum z_i_j = d_i
					
					linear = new Linear();
					linear.add(esr.getMinCapacity(),"e_"+i+"_"+j);
					linear.add(-1,"z_"+i+"_"+j);
					problem.add(new Constraint("", linear, "<=",0));
					
					linear = new Linear();
					linear.add(esr.getMaxCapacity(),"e_"+i+"_"+j);
					linear.add(-1,"z_"+i+"_"+j);
					problem.add(new Constraint("", linear, ">=",0));
										
					linear2.add(1,"e_"+i+"_"+j);
					linear3.add(1,"z_"+i+"_"+j);
					
					j++;
				}
				
				problem.add(new Constraint("", linear2, "=",1));
				
				linear3.add(-1,"d_"+i);
				problem.add(new Constraint("", linear3, "=",0));
			}
		}	
	}
	
	//add constraints related to the ESR structure in the network
	//ensures ESR configuration is valid
	//includes constraints related to selection of cost row table for each ESR
	//remove s_i_j variables
	private void setEsrConstraints_gen10() throws Exception{		
		Linear linear,linear2;
						
		for(Node node : nodes.values()){
			int i = node.getNodeID();
			
			if(node.getAllowESR()){
				double demand = node.getRequiredCapacity(esrCapacityFactor);
				for(Node n : node.getDownstreamNodes()){
					demand += n.getRequiredCapacity(esrCapacityFactor); 
				}
				int l = 0;
				for(EsrCost e: esrCost){
					if(e.getMinCapacity() > demand){
						linear = new Linear();
						linear.add(1, "e_"+i+"_"+l);
						problem.add(new Constraint("", linear, "=",0));
					}
					l++;
				}
			}
									
			linear2 = new Linear();
			linear2.add(node.getRequiredCapacity(esrCapacityFactor), "besr_"+i);
						
			for(Pipe pipe : node.getOutgoingPipes()){
				
				int j = pipe.getPipeID();
								
				//d_i sum
				double D_j = 0;
				D_j = pipe.getEndNode().getRequiredCapacity(esrCapacityFactor);
				for(Node n: pipe.getEndNode().getDownstreamNodes()){
					D_j += n.getRequiredCapacity(esrCapacityFactor);
				}
				linear2.add(D_j, "k_"+i+"_"+j);
				
				//following constraint represents besr_i = 0 => f_j = 0 where j is outgoing pipe of i 
				linear = new Linear();
				linear.add(1, "besr_"+i); 
				linear.add(-1, "f_"+j);
				problem.add(new Constraint("", linear, ">=",0));
			}
			linear2.add(-1,"d_"+i);
			problem.add(new Constraint("", linear2, "=",0));
						
			if(node.getAllowESR()){
				int j=0;
				linear2 = new Linear();
				Linear linear3 = new Linear();
				for(EsrCost esr : esrCost){
					
					//to determine e_i_j and z_i_j
					// e_i_j*min_j <= z_i_j <= e_i_j*max_j
					// sum e_i_j = 1
					// sum z_i_j = d_i
					
					linear = new Linear();
					linear.add(esr.getMinCapacity(),"e_"+i+"_"+j);
					linear.add(-1,"z_"+i+"_"+j);
					problem.add(new Constraint("", linear, "<=",0));
					
					linear = new Linear();
					linear.add(esr.getMaxCapacity(),"e_"+i+"_"+j);
					linear.add(-1,"z_"+i+"_"+j);
					problem.add(new Constraint("", linear, ">=",0));
										
					linear2.add(1,"e_"+i+"_"+j);
					linear3.add(1,"z_"+i+"_"+j);
					
					j++;
				}
				
				problem.add(new Constraint("", linear2, "=",1));
				
				linear3.add(-1,"d_"+i);
				problem.add(new Constraint("", linear3, "=",0));
			}
		}	
	}
	
	//add constraints related to pump
	private void setPumpConstraints() throws Exception{
		if(pumpGeneralProperties.pump_enabled){
			Linear linear;
			for(Pipe pipe : pipes.values()){
				int i = pipe.getPipeID();
				//f_pumphead_i = f_i * pumphead_i
				linear = new Linear();
				linear.add(1, "f_pumphead_"+i);
				linear.add(-1*maxPumpHead, "f_"+i);
				problem.add(new Constraint("", linear, "<=", 0));
				
				linear = new Linear();
				linear.add(1, "f_pumphead_"+i);
				linear.add(-1, "pumphead_"+i);
				problem.add(new Constraint("", linear, "<=", 0));
				
				linear = new Linear();
				linear.add(1, "f_pumphead_"+i);
				linear.add(-1, "pumphead_"+i);
				linear.add(-1*maxPumpHead, "f_"+i);
				problem.add(new Constraint("", linear, ">=", -1*maxPumpHead));
				
				// flow_m3h = flow*3.6;
				//power = density of water (1000) * g (9.80665) * flow (in m3h = flow*3.6) / (3.6*10^6) * efficiency
				// power = 9.81 * flow (in lps) / (1000 * efficiency)
				
				double primaryPowerCoefficient = 9.80665*pipe.getFlow()/(10*pumpGeneralProperties.efficiency);
				double secondaryPowerCoefficient = primaryPowerCoefficient*secondaryFlowFactor;
				
				// power = primarypower + secondarypower
				// power = ppc*f_i*pumphead_i + spc*(1-f_i)*pumphead_i
				// primarypower = ppc*f_pumphead_i
				// secondarypower = spc*pumphead_i - spc*f_pumphead_i
				
				linear = new Linear();
				linear.add(primaryPowerCoefficient,"f_pumphead_"+i);
				linear.add(-1,"ppumppower_"+i);
				problem.add(new Constraint("", linear, "=", 0));
				
				linear = new Linear();
				linear.add(secondaryPowerCoefficient,"pumphead_"+i);
				linear.add(-1*secondaryPowerCoefficient,"f_pumphead_"+i);
				linear.add(-1,"spumppower_"+i);
				problem.add(new Constraint("", linear, "=", 0));
				
				linear = new Linear();
				linear.add(1,"ppumppower_"+i);
				linear.add(1,"spumppower_"+i);
				linear.add(-1,"pumppower_"+i);
				problem.add(new Constraint("", linear, "=", 0));
				
				linear = new Linear();
				linear.add(1,"pumppower_"+i);
				linear.add(-1*minPumpPower,"pumphelper_"+i);
				problem.add(new Constraint("", linear, ">=", 0));
				
				linear = new Linear();
				linear.add(1,"pumppower_"+i);
				linear.add(-1*maxPumpPower,"pumphelper_"+i);
				problem.add(new Constraint("", linear, "<=", 0));
			}
		}
	}
	
	//add constraints related to pump for pipe only optimization
	private void setPumpConstraints_gen0() throws Exception{
		if(pumpGeneralProperties.pump_enabled){
			Linear linear;
			for(Pipe pipe : pipes.values()){
				int i = pipe.getPipeID();
				
				double primaryPowerCoefficient = 9.80665*pipe.getFlow()/(10*pumpGeneralProperties.efficiency);
								
				linear = new Linear();
				linear.add(primaryPowerCoefficient,"pumphead_"+i);
				linear.add(-1,"pumppower_"+i);
				problem.add(new Constraint("", linear, "=", 0));
								
				linear = new Linear();
				linear.add(1,"pumppower_"+i);
				linear.add(-1*minPumpPower,"pumphelper_"+i);
				problem.add(new Constraint("", linear, ">=", 0));
				
				linear = new Linear();
				linear.add(1,"pumppower_"+i);
				linear.add(-1*maxPumpPower,"pumphelper_"+i);
				problem.add(new Constraint("", linear, "<=", 0));
			}
		}
	}
	
	//set constraints relating ESR to pipe network
	private void setEsrPipeConstraints() throws Exception{
				
		Linear linear;
		
		for(Pipe pipe : pipes.values()){
			int i = pipe.getPipeID();
			Node end = pipe.getEndNode();
			
			//following constraint captures :
			// if node j has demand f_i = s_j_j 
			// else f_i = f_k where k is downstreampipe 
			linear = new Linear();
			linear.add(1,"f_"+i);
			
			if(end.getDemand()>0){
				int j = end.getNodeID();
				linear.add(-1,"s_"+j+"_"+j);
			}
			else{
				List<Pipe> childPipes = end.getOutgoingPipes();
				if(childPipes.size()>0){
					int k = childPipes.get(0).getPipeID();
					linear.add(-1,"f_"+k);
				}
				else
					throw new Exception("ERROR: Terminal node has no demand.");
			}
			problem.add(new Constraint("", linear, "=",0));			
		}
		
		
		for(Pipe pipe : pipes.values()){
			int i = pipe.getPipeID();
			
			for(int j=0; j<pipeCost.size();j++){
				
				//to capture y_i_j = f_i * l_i_j
				linear = new Linear();
				linear.add(pipe.getLength(),"f_"+i);
				linear.add(-1,"y_"+i+"_"+j);
				problem.add(new Constraint("", linear, ">=",0));
				
				linear = new Linear();
				linear.add(1,"l_"+i+"_"+j);
				linear.add(-1,"y_"+i+"_"+j);
				problem.add(new Constraint("", linear, ">=",0));
				
				linear = new Linear();
				linear.add(pipe.getLength(),"f_"+i);
				linear.add(-1,"y_"+i+"_"+j);
				linear.add(1,"l_"+i+"_"+j);
				problem.add(new Constraint("", linear, "<=",pipe.getLength()));
				
				if(pipe.isAllowParallel()){
					//to capture yp_i_j = f_i * p_i_j ; all boolean
					linear = new Linear();
					linear.add(1,"f_"+i);
					linear.add(-1,"yp_"+i+"_"+j);
					problem.add(new Constraint("", linear, ">=",0));
					
					linear = new Linear();
					linear.add(1,"p_"+i+"_"+j);
					linear.add(-1,"yp_"+i+"_"+j);
					problem.add(new Constraint("", linear, ">=",0));
					
					linear = new Linear();
					linear.add(1,"p_"+i+"_"+j);
					linear.add(1,"f_"+i);
					linear.add(-1,"yp_"+i+"_"+j);
					problem.add(new Constraint("", linear, "<=",1));
				}
			}
			
			if(pipe.isAllowParallel()){
				//to capture yp_i = f_i * p_i ; all boolean
				linear = new Linear();
				linear.add(1,"f_"+i);
				linear.add(-1,"yp_"+i);
				problem.add(new Constraint("", linear, ">=",0));
				
				linear = new Linear();
				linear.add(1,"p_"+i);
				linear.add(-1,"yp_"+i);
				problem.add(new Constraint("", linear, ">=",0));
				
				linear = new Linear();
				linear.add(1,"p_"+i);
				linear.add(1,"f_"+i);
				linear.add(-1,"yp_"+i);
				problem.add(new Constraint("", linear, "<=",1));
			}
		}
		
		//for all 0 demand nodes all child pipes should have same flow type f_i
		for(Node node : nodes.values()){
			if(node.getDemand()==0){
				List<Pipe> childPipes = node.getOutgoingPipes();
				if(childPipes.size()>0){
					int i = childPipes.get(0).getPipeID();
					for(Pipe pipe : childPipes){
						int j = pipe.getPipeID();
						if(i!=j){
							linear = new Linear();
							linear.add(1,"f_"+i);
							linear.add(-1,"f_"+j);
							problem.add(new Constraint("", linear, "=",0));
							i=j;
						}
					}
				}
				else
					throw new Exception("ERROR: Terminal node has no demand");
			}
		}
		
		//if s_i_j=1 then all pipes in path from i to j should all be 0
		for(Node node : nodes.values()){
			if(node.getDemand()>0){
				List<Pipe> sourceToNodePipes = node.getSourceToNodePipes();
				int j = node.getNodeID();
				for(int i=0;i<sourceToNodePipes.size();i++){
					Node startNode = sourceToNodePipes.get(i).getStartNode();
					if(startNode.getDemand()>0){
						int sid = startNode.getNodeID();
						for(int k=i;k<sourceToNodePipes.size();k++){
							int id = sourceToNodePipes.get(k).getPipeID();
							linear = new Linear();
							linear.add(1,"s_"+sid+"_"+j);
							linear.add(1,"f_"+id);
							problem.add(new Constraint("", linear, "<=",1));
						}
					}
				}
			}
		}
		
//		//if s_i_j=1 then first pipe in path from i to j should be 0
//		for(Node node : nodes.values()){
//			if(node.getDemand()>0){
//				List<Pipe> sourceToNodePipes = node.getSourceToNodePipes();
//				int j = node.getNodeID();
//				for(int i=0;i<sourceToNodePipes.size();i++){
//					Pipe pipe = sourceToNodePipes.get(i);
//					Node startNode = pipe.getStartNode();
//					if(startNode.getDemand()>0){
//						int sid = startNode.getNodeID();
//						
//						int id = pipe.getPipeID();
//						linear = new Linear();
//						linear.add(1,"s_"+sid+"_"+j);
//						linear.add(1,"f_"+id);
//						problem.add(new Constraint("", linear, "<=",1));
//
//					}
//				}
//			}
//		}
		
	}
	
	//set constraints relating ESR to pipe network
	//allow zero demand nodes to have ESRs
	private void setEsrPipeConstraints_gen() throws Exception{
		
		Linear linear;
		
		for(Pipe pipe : pipes.values()){
			int i = pipe.getPipeID();
			Node end = pipe.getEndNode();
			
			//following constraint captures :
			// f_i = s_j_j 
			linear = new Linear();
			linear.add(1,"f_"+i);
			
			int j = end.getNodeID();
			linear.add(-1,"s_"+j+"_"+j);
			
			problem.add(new Constraint("", linear, "=",0));			
		}
		
		
		for(Pipe pipe : pipes.values()){
			int i = pipe.getPipeID();
			
			for(int j=0; j<pipeCost.size();j++){
				
				//to capture y_i_j = f_i * l_i_j
				linear = new Linear();
				linear.add(pipe.getLength(),"f_"+i);
				linear.add(-1,"y_"+i+"_"+j);
				problem.add(new Constraint("", linear, ">=",0));
				
				linear = new Linear();
				linear.add(1,"l_"+i+"_"+j);
				linear.add(-1,"y_"+i+"_"+j);
				problem.add(new Constraint("", linear, ">=",0));
				
				linear = new Linear();
				linear.add(pipe.getLength(),"f_"+i);
				linear.add(-1,"y_"+i+"_"+j);
				linear.add(1,"l_"+i+"_"+j);
				problem.add(new Constraint("", linear, "<=",pipe.getLength()));
				
				if(pipe.isAllowParallel()){
					//to capture yp_i_j = f_i * p_i_j ; all boolean
					linear = new Linear();
					linear.add(1,"f_"+i);
					linear.add(-1,"yp_"+i+"_"+j);
					problem.add(new Constraint("", linear, ">=",0));
					
					linear = new Linear();
					linear.add(1,"p_"+i+"_"+j);
					linear.add(-1,"yp_"+i+"_"+j);
					problem.add(new Constraint("", linear, ">=",0));
					
					linear = new Linear();
					linear.add(1,"p_"+i+"_"+j);
					linear.add(1,"f_"+i);
					linear.add(-1,"yp_"+i+"_"+j);
					problem.add(new Constraint("", linear, "<=",1));
				}
			}
			
			if(pipe.isAllowParallel()){
				//to capture yp_i = f_i * p_i ; all boolean
				linear = new Linear();
				linear.add(1,"f_"+i);
				linear.add(-1,"yp_"+i);
				problem.add(new Constraint("", linear, ">=",0));
				
				linear = new Linear();
				linear.add(1,"p_"+i);
				linear.add(-1,"yp_"+i);
				problem.add(new Constraint("", linear, ">=",0));
				
				linear = new Linear();
				linear.add(1,"p_"+i);
				linear.add(1,"f_"+i);
				linear.add(-1,"yp_"+i);
				problem.add(new Constraint("", linear, "<=",1));
			}
		}
	}
	
	//set constraints relating ESR to pipe network
	//allow zero demand nodes to have ESRs
	//some constraints removed due to introduction of l_i_j_k
	private void setEsrPipeConstraints_gen3() throws Exception{
		
		Linear linear;
		
		for(Pipe pipe : pipes.values()){
			int i = pipe.getPipeID();
			Node end = pipe.getEndNode();
			
			//following constraint captures :
			// f_i = s_j_j 
			linear = new Linear();
			linear.add(1,"f_"+i);
			
			int j = end.getNodeID();
			linear.add(-1,"s_"+j+"_"+j);
			
			problem.add(new Constraint("", linear, "=",0));			
		}		
	}
	
	//set constraints relating ESR to pipe network
	//allow zero demand nodes to have ESRs
	//some constraints removed due to introduction of l_i_j_k
	//remove s_i_j variables
	private void setEsrPipeConstraints_gen4() throws Exception{
		
		Linear linear;
		
		for(Pipe pipe : pipes.values()){
			int i = pipe.getPipeID();
			Node end = pipe.getEndNode();
			
			//following constraint captures :
			// f_i = s_j_j 
			linear = new Linear();
			linear.add(1,"f_"+i);
			
			int j = end.getNodeID();
			linear.add(-1,"besr_"+j);
			
			problem.add(new Constraint("", linear, "=",0));			
		}		
		
		for(Pipe pipe: source.getOutgoingPipes()){
			int i = pipe.getPipeID();
			linear = new Linear();
			linear.add(1,"f_"+i);
			
			problem.add(new Constraint("", linear, "=",1));
		}
	}
	

	//post optimization, set the water head at different nodes
	private void setHeadsFromResult(Result result) throws Exception{
		int j=0;
		for(Node node : nodes.values()){
			double headloss = 0;
			for(Pipe pipe : node.getSourceToNodePipes()){
				int i = pipe.getPipeID();
				if(pipe.isAllowParallel()){
					j=0;
					//double temp = 0;
					for(PipeCost entry : pipeCost){				
						// primary flow in case of parallel pipe
						double flow = pipe.getFlow() / (1 + (entry.getRoughness()/pipe.getRoughness())*Math.pow(entry.getDiameter()/pipe.getDiameter(), 4.87/1.852));
						double loss = Util.HWheadLoss(pipe.getLength(), flow, pipe.getRoughness(), pipe.getDiameter());
						//double temp2 = result.getPrimalValue(("p_"+i+"_"+j)).doubleValue();
						//temp += temp2;
						headloss += loss * result.getPrimalValue(("p_"+i+"_"+j)).intValue();
						j++;
					}	
					//double temp2 = result.getPrimalValue(("p_"+i)).doubleValue();
					//temp += temp2;
					//System.out.println(pipe.getPipeID() + " " + temp);
					double flow = pipe.getFlow();
					double loss = Util.HWheadLoss(pipe.getLength(), flow, pipe.getRoughness(), pipe.getDiameter());
					headloss += loss * result.getPrimalValue(("p_"+i)).intValue();
				}
				else{
					if(pipe.getDiameter()!=0){
						j=0;
						for(PipeCost entry : pipeCost){
							double loss = Util.HWheadLoss(pipe.getFlow(), pipe.getRoughness(), entry.getDiameter()); 
							double length = result.getPrimalValue(("l_"+i+"_"+j)).doubleValue();
							length = Util.round(length, 5);
							headloss += loss * length;
							j++;
						}
					}
					else{
						j=0;
						for(PipeCost entry : pipeCost){
							double loss = Util.HWheadLoss(pipe.getFlow(), entry.getRoughness(), entry.getDiameter()); 
							double length = result.getPrimalValue(("l_"+i+"_"+j)).doubleValue();
							length = Util.round(length, 5);
							headloss += loss * length;
							j++;
						}	
					}
				}
				
				if(pumpGeneralProperties.pump_enabled){
					double pumphead = result.getPrimalValue("pumphead_"+i).doubleValue();
					headloss = headloss - pumphead;
				}
				
				headloss = headloss + pipe.getValveSetting();
			}
			//ESR height
//			if(node.getDemand()>0){
//				double esr_height = result.getPrimalValue(("h_"+node.getNodeID())).doubleValue();
//				System.out.println(("h_"+node.getNodeID())+" "+esr_height);
//			}
			node.setHead(source.getHead() - headloss);
		}
	}
			

	//get the ESR cost, given the capacity in litres
	public double getCost(double capacity, List<EsrCost> esrCosts){
		for(EsrCost e: esrCosts){
			if(e.getMinCapacity() <= capacity && e.getMaxCapacity() > capacity){
				return e.getBaseCost() + (capacity-e.getMinCapacity())*e.getUnitCost();
			}
		}
		return Integer.MAX_VALUE;
	}
	
	//after optimization set water head for nodes 
	private void setHeadsFromResult_esr(Result result) throws Exception{

		for(Node node : nodes.values()){
			double head = result.getPrimalValue("head_"+node.getNodeID()).doubleValue();
			node.setHead(head);
			int i = node.getNodeID();
			int esrChoice = result.getPrimalValue("s_"+i+"_"+i).intValue();
			double demandServed = result.getPrimalValue("d_"+i).doubleValue();
			demandServed = Util.round(demandServed, 5);
			if(esrChoice==1 && demandServed > 0){
				node.setESR(i);
				double esrHeight = result.getPrimalValue("esr_"+i).doubleValue();
				node.setEsrHeight(esrHeight);
				node.setEsrTotalDemand(demandServed);
				node.setEsrCost(getCost(demandServed, esrCost));
				node.addToServedNodes(node);
				
				for(Node downstreamNode : node.getDownstreamNodes()){
					int j = downstreamNode.getNodeID();
					esrChoice = result.getPrimalValue("s_"+i+"_"+j).intValue();
					if(esrChoice==1){
						downstreamNode.setESR(i);
						node.addToServedNodes(downstreamNode);
					}
				}
			}
		}
	}
	
	//after optimization set water head for nodes
	//allow zero demand nodes to have ESRs 
	private void setHeadsFromResult_esr_gen2(Result result) throws Exception{

		for(Node node : nodes.values()){
			double head = result.getPrimalValue("head_"+node.getNodeID()).doubleValue();
			node.setHead(head);
			int i = node.getNodeID();
			int esrChoice = result.getPrimalValue("besr_"+i).intValue();
			double demandServed = result.getPrimalValue("d_"+i).doubleValue();
			demandServed = Util.round(demandServed, 5);
			if(esrChoice==1 && demandServed > 0){
				node.setESR(i);
				double esrHeight = result.getPrimalValue("esr_"+i).doubleValue();
				node.setEsrHeight(esrHeight);
				node.setEsrTotalDemand(demandServed);
				node.setEsrCost(getCost(demandServed, esrCost));
				node.addToServedNodes(node);
				
				for(Pipe pipe : node.getOutgoingPipes()){	
					int j = pipe.getPipeID();
					esrChoice = 1 - result.getPrimalValue("f_"+j).intValue();
					if(esrChoice==1){
						Node endnode = pipe.getEndNode();
						endnode.setESR(i);
						node.addToServedNodes(endnode);
						
						for(Node n : endnode.getDownstreamNodes()){
							n.setESR(i);
							node.addToServedNodes(n);
						}
					}
				}
			}
		}
	}
	

		
	//after optimization set diameters for pipes
	private void setDiametersFromResult(Result result) throws Exception {
		int j;
		for(Pipe pipe : pipes.values()){
			j=0;
			int i = pipe.getPipeID();
			int noOfSubPipes = 0;
			
			if(pumpGeneralProperties.pump_enabled){
				double pumphead = result.getPrimalValue("pumphead_"+i).doubleValue();
				double power = result.getPrimalValue("pumppower_"+i).doubleValue();
				
				pumphead = Util.round(pumphead, 5);
				power = Util.round(power, 5);
				
				pipe.setPumpHead(pumphead);
				pipe.setPumpPower(power);
				
				double presentvaluefactor = Util.presentValueFactor(pumpGeneralProperties.discount_rate, pumpGeneralProperties.inflation_rate, pumpGeneralProperties.design_lifetime);		
				double primarycoeffecient = 365*presentvaluefactor*generalProperties.supply_hours*pumpGeneralProperties.energycost_per_kwh*pumpGeneralProperties.energycost_factor;
				double energycost = power*primarycoeffecient;
				double capitalcost = power*pumpGeneralProperties.capitalcost_per_kw;
				
				System.out.println("pipe:"+i+" pumphead:"+pumphead);
				System.out.println("Power: "+power+" Energy Cost: "+energycost+" Capital Cost:" + capitalcost);				
			}
			
			
			if(pipe.existingPipe()){
				for(PipeCost entry : pipeCost){
					if(pipe.isAllowParallel()){
						double parallelChoice = result.getPrimalValue(("p_"+i+"_"+j)).intValue();
						if(parallelChoice==1){
							pipe.setDiameter2(entry.getDiameter());
							pipe.setRoughness2(entry.getRoughness());
							pipe.setChosenPipeCost2(entry);
						}
					}		
					double length = result.getPrimalValue(("l_"+i+"_"+j)).doubleValue();
					length = Util.round(length, 5);
					if(entry.getDiameter() == pipe.getDiameter()){	
						if(length != pipe.getLength())
							throw new Exception("Something wrong in parallel link "+i);
					}
					j++;
				}
			}
			else{
				for(PipeCost entry : pipeCost){
					double length = result.getPrimalValue(("l_"+i+"_"+j)).doubleValue();
					length = Util.round(length, 5);
					if(length>0){
						noOfSubPipes++;
						if(noOfSubPipes==1){
							pipe.setDiameter(entry.getDiameter());
							pipe.setRoughness(entry.getRoughness());
							pipe.setChosenPipeCost(entry);
						}
						else if(noOfSubPipes==2){
							pipe.setDiameter2(entry.getDiameter());
							pipe.setLength2(length);
							pipe.setRoughness2(entry.getRoughness());
							pipe.setChosenPipeCost2(entry);
						}
						else
							throw new Exception("Could not solve with given constraints");
					}
					j++;
				}
			}
		}
	}
	
	//after optimization set diameters for pipes
	//updated with change of l_i_j to l_i_j_k
	private void setDiametersFromResult_gen3(Result result) throws Exception{
		int j;
		for(Pipe pipe : pipes.values()){
			j=0;
			int i = pipe.getPipeID();
			int noOfSubPipes = 0;
			
			int flowChoice = result.getPrimalValue("f_"+i).intValue();
			pipe.setFlowchoice(flowChoice==1?FlowType.PRIMARY:FlowType.SECONDARY);			
			
			if(pipe.existingPipe()){
				for(PipeCost entry : pipeCost){
					if(pipe.isAllowParallel()){
						double parallelChoice1 = result.getPrimalValue(("p_"+i+"_"+j+"_0")).intValue();
						double parallelChoice2 = result.getPrimalValue(("p_"+i+"_"+j+"_1")).intValue();
						
						if(parallelChoice1==1 || parallelChoice2==1){
							pipe.setDiameter2(entry.getDiameter());
							pipe.setRoughness2(entry.getRoughness());
							pipe.setChosenPipeCost2(entry);
						}
					}		
					double length1 = result.getPrimalValue(("l_"+i+"_"+j+"_0")).doubleValue();
					double length2 = result.getPrimalValue(("l_"+i+"_"+j+"_1")).doubleValue();
					
					length1 = Util.round(length1, 5);
					length2 = Util.round(length2, 5);
					double length = Math.max(length1, length2);
					if(entry.getDiameter() == pipe.getDiameter()){	
						if(length != pipe.getLength())
							throw new Exception("Something wrong in parallel link "+i);
					}
					j++;
				}
			}
			else{
				for(PipeCost entry : pipeCost){
					double length1 = result.getPrimalValue(("l_"+i+"_"+j+"_0")).doubleValue();
					double length2 = result.getPrimalValue(("l_"+i+"_"+j+"_1")).doubleValue();
					
					length1 = Util.round(length1, 5);
					length2 = Util.round(length2, 5);
					double length = Math.max(length1, length2);
					if(length>0){
						noOfSubPipes++;
						if(noOfSubPipes==1){
							pipe.setDiameter(entry.getDiameter());
							pipe.setRoughness(entry.getRoughness());
							pipe.setChosenPipeCost(entry);
						}
						else if(noOfSubPipes==2){
							pipe.setDiameter2(entry.getDiameter());
							pipe.setLength2(length);
							pipe.setRoughness2(entry.getRoughness());
							pipe.setChosenPipeCost2(entry);
						}
						else
							throw new Exception("more than 2 pipes for link "+i);
					}
					j++;
				}
			}
		}
	}
	
	//after optimization set diameters for pipes
	//included pumps
	private void setDiametersFromResult_gen4(Result result) throws Exception{
		int j;
		for(Pipe pipe : pipes.values()){
			j=0;
			int i = pipe.getPipeID();
			int noOfSubPipes = 0;
			
			int flowChoice = result.getPrimalValue("f_"+i).intValue();
			pipe.setFlowchoice(flowChoice==1?FlowType.PRIMARY:FlowType.SECONDARY);			
			
			if(pumpGeneralProperties.pump_enabled){
				double pumphead = result.getPrimalValue("pumphead_"+i).doubleValue();
				double primarypower = result.getPrimalValue("ppumppower_"+i).doubleValue();
				double secondarypower = result.getPrimalValue("spumppower_"+i).doubleValue();
				double power = result.getPrimalValue("pumppower_"+i).doubleValue();
				
				pumphead = Util.round(pumphead, 5);
				primarypower = Util.round(primarypower, 5);
				secondarypower = Util.round(secondarypower, 5);
				power = Util.round(power, 5);
				
				pipe.setPumpHead(pumphead);
				pipe.setPumpPower(power);
				
				double presentvaluefactor = Util.presentValueFactor(pumpGeneralProperties.discount_rate, pumpGeneralProperties.inflation_rate, pumpGeneralProperties.design_lifetime);		
				double primarycoeffecient = 365*presentvaluefactor*generalProperties.supply_hours*pumpGeneralProperties.energycost_per_kwh*pumpGeneralProperties.energycost_factor;
				double secondarycoeffecient = 365*presentvaluefactor*esrGeneralProperties.secondary_supply_hours*pumpGeneralProperties.energycost_per_kwh*pumpGeneralProperties.energycost_factor;
				double energycost = primarypower*primarycoeffecient + secondarypower*secondarycoeffecient;
				double capitalcost = power*pumpGeneralProperties.capitalcost_per_kw;
				
				System.out.println("pipe:"+i+" pumphead:"+pumphead);
				System.out.println("Power: "+power+" Energy Cost: "+energycost+" Capital Cost:" + capitalcost);				
			}
			
			if(pipe.existingPipe()){
				for(PipeCost entry : pipeCost){
					if(pipe.isAllowParallel()){
						double parallelChoice1 = result.getPrimalValue(("p_"+i+"_"+j+"_0")).intValue();
						double parallelChoice2 = result.getPrimalValue(("p_"+i+"_"+j+"_1")).intValue();
						
						if(parallelChoice1==1 || parallelChoice2==1){
							pipe.setDiameter2(entry.getDiameter());
							pipe.setRoughness2(entry.getRoughness());
							pipe.setChosenPipeCost2(entry);
						}
					}		
					double length1 = result.getPrimalValue(("l_"+i+"_"+j+"_0")).doubleValue();
					double length2 = result.getPrimalValue(("l_"+i+"_"+j+"_1")).doubleValue();
					
					length1 = Util.round(length1, 5);
					length2 = Util.round(length2, 5);
					double length = Math.max(length1, length2);
					if(entry.getDiameter() == pipe.getDiameter()){	
						if(length != pipe.getLength())
							throw new Exception("Something wrong in parallel link "+i);
					}
					j++;
				}
			}
			else{
				for(PipeCost entry : pipeCost){
					double length1 = result.getPrimalValue(("l_"+i+"_"+j+"_0")).doubleValue();
					double length2 = result.getPrimalValue(("l_"+i+"_"+j+"_1")).doubleValue();
					
					length1 = Util.round(length1, 5);
					length2 = Util.round(length2, 5);
					double length = Math.max(length1, length2);
					if(length>0){
						noOfSubPipes++;
						if(noOfSubPipes==1){
							pipe.setDiameter(entry.getDiameter());
							pipe.setRoughness(entry.getRoughness());
							pipe.setChosenPipeCost(entry);
						}
						else if(noOfSubPipes==2){
							pipe.setDiameter2(entry.getDiameter());
							pipe.setLength2(length);
							pipe.setRoughness2(entry.getRoughness());
							pipe.setChosenPipeCost2(entry);
						}
						else
							throw new Exception("more than 2 pipes for link "+i);
					}
					j++;
				}
			}
		}
	}
	



	// ---
	public static int attempts=0;
	// This should point to the root of the repository: Jaltantra-Code-and-Scripts
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

	static String PREVIOUS_SOLVER_EXECUTION_TIME=SOLVER_EXECUTION_TIME;
	static String SOLVER_EXECUTION_TIME_DISPLAY_STR = "5 minutes";

	static String list_of_times[]={"00:10:00","00:20:00"};

	/**
	 * Optimize the network. And, if done successfully, the results are stored in three ArrayLists,
	 * namely resultPipes, resultCost and resultPumps.
	 *
	 * @return whether network was solved successfully results are ready (`true`) or not (`false`)
	 * @throws Exception in case of any error/problem or to convey any status information
	 */
	public boolean Optimize() throws Exception {
		// Execution flow:
		//   1. Create the data files for the network
		//   2. Asynchronously launch `CalculateNetworkCost.py` for the data files

		// Validate the network layout
		logd("Network validation started...");
		final int networkValidationResult = validateNetwork();
		logd("Network validation complete..., networkValidationResult = " + networkValidationResult);
		if (networkValidationResult == 1 || networkValidationResult == 2) {
			final String networkFileResult = createNetworkFile();
			final String networkFileStatus = networkFileResult.substring(0, networkFileResult.indexOf("-"));
			final String networkFileName = networkFileResult.substring(2);
			final String networkFileHash = networkFileName.substring(0, networkFileName.lastIndexOf("."));
			boolean statusFileExists = (new File(SOLVER_ROOT_DIR + "/" + SOLVER_3_AUTO_SOLVE_SCRIPT_DIR + "/" + networkFileHash + "/0_status")).exists();
			logd("networkFileHash = " + networkFileHash + ", statusFileExists = " + statusFileExists);

			if (networkFileStatus.equals("0") || statusFileExists == false) {
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
					throw new Exception("Internal server error: failed to start the solver launch script");
				} else if (status == 0) {
					loge("FIXME: Failed to start the solver for network file: " + networkFileResult);
					throw new Exception("Internal server error: failed to start the solver for this network");
				} else if (status == 1) {
					String filePath=SOLVER_ROOT_DIR+"/"+SOLVER_3_AUTO_SOLVE_SCRIPT_DIR+"/"+networkFileHash+"/baron_m1_"+networkFileHash+"/"+"std_out_err.txt";
					File file = new File(filePath);
					if(file.exists()){
						try (BufferedReader reader = new BufferedReader(new FileReader(file))) {
							String line = reader.readLine();
							System.out.println("File content: " + line);
						} catch (IOException e) {
							System.out.println("Error reading from file: " + e.getMessage());
						}
						System.out.println("File readed successfully.");
					}
					logi("Solver is running for this network. Please wait...");
					throw new Exception("Solver is running for this network. Please wait...");
				} else if (status == 2) {
					loge("CHECKME: The solvers finished the execution, but failed to get the result. " +
						 "Either some unknown error, or no feasible solution found, or failed to solve the network, click optimize again and check after " + SOLVER_EXECUTION_TIME_DISPLAY_STR+" minutes");
					//creating an attempts file
					String pathToattemptsFile=SOLVER_ROOT_DIR+"/"+SOLVER_2_HASH_FILE_DIR+"/"+networkFileHash;
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
						System.out.println("system is now executing for 10 minutes");
						if (status == 3) {
							logi("Extracting the result");
							boolean ok = extractSolverResult(SOLVER_ROOT_DIR + "/" + SOLVER_3_AUTO_SOLVE_SCRIPT_DIR + "/" + networkFileHash);
							if (ok) return true;
							loge("CHECKME: FIXME: extractSolverResult(...) return false for network file with hash =" + networkFileHash);
							throw new Exception("Internal server error: result extraction failed for the network file with hash =" + networkFileHash);
						}
						SOLVER_EXECUTION_TIME_DISPLAY_STR=""+getMinutesToDisplay(SOLVER_EXECUTION_TIME);
						throw new Exception("Either no feasible solution found, or failed to solve the network, retrying for " + SOLVER_EXECUTION_TIME_DISPLAY_STR+" minutes");
					}
					throw new Exception("Either no feasible solution found, or failed to solve the network in " + SOLVER_EXECUTION_TIME_DISPLAY_STR+" minutes");

					// if(attempts == 1){
					// 	attempts=2;
					// 	try (FileWriter writer = new FileWriter(file)) {
					// 		writer.write(""+attempts);
					// 		System.out.println("File written successfully.");
					// 	} catch (IOException e) {
					// 		System.out.println("Error writing to file: " + e.getMessage());
					// 	}
					// 	SOLVER_EXECUTION_TIME="00:10:00";
					// 	System.out.println("system is running for "+SOLVER_EXECUTION_TIME+"minutes");
					//
					// 	//delete previously existed network results
					// 	String pathToFile=SOLVER_ROOT_DIR + "/" + SOLVER_3_AUTO_SOLVE_SCRIPT_DIR + "/" + networkFileHash;
					// 	deleteNetworkResults(pathToFile);
					// 	statusFileExists = (new File(SOLVER_ROOT_DIR + "/" + SOLVER_3_AUTO_SOLVE_SCRIPT_DIR + "/" + networkFileHash + "/0_status")).exists();
					// 	launchCalculateNetworkCost(SOLVER_ROOT_DIR + "/" + SOLVER_2_HASH_FILE_DIR + "/" + networkFileName);
					// 	status = checkSolverResultStatus(SOLVER_ROOT_DIR + "/" + SOLVER_3_AUTO_SOLVE_SCRIPT_DIR + "/" + networkFileHash);
					// 	System.out.println("system is now executing for 10 minutes");
					// 	if (status == 3) {
					// 		logi("Extracting the result");
					// 		boolean ok = extractSolverResult(SOLVER_ROOT_DIR + "/" + SOLVER_3_AUTO_SOLVE_SCRIPT_DIR + "/" + networkFileHash);
					// 		if (ok) return true;
					// 		loge("CHECKME: FIXME: extractSolverResult(...) return false for network file with hash =" + networkFileHash);
					// 		throw new Exception("Internal server error: result extraction failed for the network file with hash =" + networkFileHash);
					// 	}
					// 	SOLVER_EXECUTION_TIME_DISPLAY_STR=""+getMinutesToDisplay(SOLVER_EXECUTION_TIME);
					// 	throw new Exception("Either no feasible solution found, or failed to solve the network, retrying for " + SOLVER_EXECUTION_TIME_DISPLAY_STR+" minutes");
					// }else
					// if(attempts == 2){
					// 	attempts=3;
					// 	try (FileWriter writer = new FileWriter(file)) {
					// 		writer.write(""+attempts);
					// 		System.out.println("File written successfully.");
					// 	} catch (IOException e) {
					// 		System.out.println("Error writing to file: " + e.getMessage());
					// 	}
					// 	SOLVER_EXECUTION_TIME="00:20:00";
					// 	System.out.println("system is running for "+SOLVER_EXECUTION_TIME+"minutes");
					//
					// 	//delete previously existed network results
					// 	String pathToFile=SOLVER_ROOT_DIR + "/" + SOLVER_3_AUTO_SOLVE_SCRIPT_DIR + "/" + networkFileHash;
					// 	deleteNetworkResults(pathToFile);
					// 	statusFileExists = (new File(SOLVER_ROOT_DIR + "/" + SOLVER_3_AUTO_SOLVE_SCRIPT_DIR + "/" + networkFileHash + "/0_status")).exists();
					// 	launchCalculateNetworkCost(SOLVER_ROOT_DIR + "/" + SOLVER_2_HASH_FILE_DIR + "/" + networkFileName);
					// 	status = checkSolverResultStatus(SOLVER_ROOT_DIR + "/" + SOLVER_3_AUTO_SOLVE_SCRIPT_DIR + "/" + networkFileHash);
					// 	System.out.println("system is running for 20 minutes");
					// 	if (status == 3) {
					// 		logi("Extracting the result");
					// 		boolean ok = extractSolverResult(SOLVER_ROOT_DIR + "/" + SOLVER_3_AUTO_SOLVE_SCRIPT_DIR + "/" + networkFileHash);
					// 		if (ok) return true;
					// 		loge("CHECKME: FIXME: extractSolverResult(...) return false for network file with hash =" + networkFileHash);
					// 		throw new Exception("Internal server error: result extraction failed for the network file with hash =" + networkFileHash);
					// 	}
					// 	SOLVER_EXECUTION_TIME_DISPLAY_STR=""+getMinutesToDisplay(SOLVER_EXECUTION_TIME);
					// 	throw new Exception("Either no feasible solution found, or failed to solve the network, retrying for " + SOLVER_EXECUTION_TIME_DISPLAY_STR+" minutes");
					// }else throw new Exception("Either no feasible solution found, or failed to solve the network, in " + SOLVER_EXECUTION_TIME_DISPLAY_STR+" minutes");

					// String timeOptions[]={"00:10:00","01:00:00","04:00:00",};
					// for(String str:timeOptions) {
					// 	SOLVER_EXECUTION_TIME = str;
					// 	System.out.println("system is running for " + SOLVER_EXECUTION_TIME + "minutes");
					//
					// 	//delete previously existed network results
					// 	String pathToFile = SOLVER_ROOT_DIR + "/" + SOLVER_3_AUTO_SOLVE_SCRIPT_DIR + "/" + networkFileHash;
					// 	deleteNetworkResults(pathToFile);
					// 	statusFileExists = (new File(SOLVER_ROOT_DIR + "/" + SOLVER_3_AUTO_SOLVE_SCRIPT_DIR + "/" + networkFileHash + "/0_status")).exists();
					// 	launchCalculateNetworkCost(SOLVER_ROOT_DIR + "/" + SOLVER_2_HASH_FILE_DIR + "/" + networkFileName);
					// 	status = checkSolverResultStatus(SOLVER_ROOT_DIR + "/" + SOLVER_3_AUTO_SOLVE_SCRIPT_DIR + "/" + networkFileHash);
					// 	System.out.println("system got executed for 20 minutes");
					// 	if (status == 3) {
					// 		logi("Extracting the result");
					// 		boolean ok = extractSolverResult(SOLVER_ROOT_DIR + "/" + SOLVER_3_AUTO_SOLVE_SCRIPT_DIR + "/" + networkFileHash);
					// 		if (ok) return true;
					// 		loge("CHECKME: FIXME: extractSolverResult(...) return false for network file with hash =" + networkFileHash);
					// 		throw new Exception("Internal server error: result extraction failed for the network file with hash =" + networkFileHash);
					// 	}
					// 	SOLVER_EXECUTION_TIME_DISPLAY_STR = "" + getMinutesToDisplay(SOLVER_EXECUTION_TIME);
					// }
					// throw new Exception("Either no feasible solution found, or failed to solve the network, in " + SOLVER_EXECUTION_TIME_DISPLAY_STR+" minutes");
				} else if (status == 3) {
					logi("Extracting the result");
					boolean ok = extractSolverResult(SOLVER_ROOT_DIR + "/" + SOLVER_3_AUTO_SOLVE_SCRIPT_DIR + "/" + networkFileHash);
					if (ok) return true;
					loge("CHECKME: FIXME: extractSolverResult(...) return false for network file with hash =" + networkFileHash);
					throw new Exception("Internal server error: result extraction failed for the network file with hash =" + networkFileHash);
				} else {
					loge("CHECKME: FIXME: Unexpected return value from `checkSolverResultStatus()`. Need to delete the past execution result for hash =" + networkFileHash);
					throw new Exception("Internal server error: unknown execution status. Need to delete the past execution result for hash =" + networkFileHash);
				}
			} else {
				loge("FIXME: `createNetworkFile()` method returned unexpected value ' " + networkFileResult + " '");
				throw new Exception("Internal server error: createNetworkFile() method failed, probably due to unexpected change in some method(s)");
			}
		} else if (networkValidationResult == 3) {
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

	public static String parseDate(String dateString,int multiplier) {
		SimpleDateFormat dateFormat = new SimpleDateFormat("HH:mm:ss");
		String updatedTime="";
		try {
			System.out.println(dateString);
			Date date = dateFormat.parse(dateString);

			int doubled=date.getMinutes() * multiplier;
			System.out.println("minutes : "+ doubled);

			int hours = doubled / 60;
			int remainingMinutes = doubled % 60;
			int seconds = 0;

			updatedTime = String.format("%02d:%02d:%02d", hours, remainingMinutes, seconds);
			System.out.println("Time: " + updatedTime);

		} catch (ParseException e) {
			System.out.println("Error parsing date string: " + e.getMessage());
		}

		return updatedTime;
	}

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
		return "0-" + hashedFileName;
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
					SOLVER_EXECUTION_TIME
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
		logi("Best result = " + resultLines[6]);

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

	private boolean acyclic_Optimize(long startTime) throws Exception
	{
		
		//switch between different models and setup downstream/upstream node information
		switch(modelNumber){
			case 0:
			case 1:
				getNodeSupply(source);
				setSourceToNodePipes(source);
				break;
			case 2:
			case 3:
			case 4:
			case 5:
			case 6:
			case 7:
			case 8:
			case 9:
			case 10:
				getNodeSupply_gen(source);
				setSourceToNodePipes_gen(source);
				break;
		}
				
		//legacy code, earlier ILP library used was GPLK, can be removed
//				SolverFactory factory = new SolverFactoryGLPK(); // use GPLK
//				
//				factory.setParameter(Solver.VERBOSE, 0); 
//				factory.setParameter(Solver.TIMEOUT, 50); // set timeout in seconds
		
		problem = new Problem();
		//problem.setTimeLimit(50000); // set timeout in milliseconds for web
		problem.setTimeLimit(60000*60*10); // set timeout in milliseconds for local
		
		//problem.setTimeLimit(4*60*60*1000);
				
		//switch between different models and set up the objective cost and constraints of the ILP model
		switch(modelNumber){
			case 0:
				setObjectiveCost();
				System.out.println("After setobjective: "+problem.getConstraintsCount());
				setPumpOptionConstraints();
				System.out.println("After setpumpoption: "+problem.getConstraintsCount());
				setValveOptionConstraints();
				System.out.println("After setvalveoption: "+problem.getConstraintsCount());
				setPipeConstraints();
				System.out.println("After pipeconstraints: "+problem.getConstraintsCount());
				setHeadLossConstraints();
				System.out.println("After headlossconstraints: "+problem.getConstraintsCount());
				setPumpConstraints_gen0();
				System.out.println("After pumpconstraints: "+problem.getConstraintsCount());
				break;
			case 1:
				addVariables();
				System.out.println("After addvariables: "+problem.getConstraintsCount());
				setObjectiveCost_esr();
				System.out.println("After setobjective: "+problem.getConstraintsCount());
				setPipeConstraints();
				System.out.println("After pipeconstraints: "+problem.getConstraintsCount());
				setHeadLossConstraints_esr();
				System.out.println("After headlossconstraints: "+problem.getConstraintsCount());
				setEsrConstraints();
				System.out.println("After esrconstraints: "+problem.getConstraintsCount());
				setEsrPipeConstraints();
				System.out.println("After esrpipeconstraints: "+problem.getConstraintsCount());
				break;
			case 2:
				addVariables_gen();
				System.out.println("After addvariables: "+problem.getConstraintsCount());
				setObjectiveCost_esr_gen();
				System.out.println("After setobjective: "+problem.getConstraintsCount());
				setPipeConstraints();
				System.out.println("After pipeconstraints: "+problem.getConstraintsCount());
				setHeadLossConstraints_esr_gen();
				System.out.println("After headlossconstraints: "+problem.getConstraintsCount());
				setEsrConstraints_gen();
				System.out.println("After esrconstraints: "+problem.getConstraintsCount());
				setEsrPipeConstraints_gen();
				System.out.println("After esrpipeconstraints: "+problem.getConstraintsCount());
				break;
			case 3:
				addVariables_gen();
				System.out.println("After addvariables:aaaa "+problem.getConstraintsCount());
				setObjectiveCost_esr_gen();
				System.out.println("After setobjective: "+problem.getConstraintsCount());
				setPipeConstraints_model3();
				System.out.println("After pipeconstraints: "+problem.getConstraintsCount());
				setHeadLossConstraints_esr_gen2();
				System.out.println("After headlossconstraints: "+problem.getConstraintsCount());
				setEsrConstraints_gen2();
				System.out.println("After esrconstraints: "+problem.getConstraintsCount());
				setEsrPipeConstraints_gen();
				System.out.println("After esrpipeconstraints: "+problem.getConstraintsCount());
				break;
			case 4:
				addVariables_gen3();
				System.out.println("After addvariables: "+problem.getConstraintsCount());
				setObjectiveCost_esr_gen3();
				System.out.println("After setobjective: "+problem.getConstraintsCount());
				setPipeConstraints_gen3();
				System.out.println("After pipeconstraints: "+problem.getConstraintsCount());
				setHeadLossConstraints_esr_gen3();
				System.out.println("After headlossconstraints: "+problem.getConstraintsCount());
				setEsrConstraints_gen2();
				System.out.println("After esrconstraints: "+problem.getConstraintsCount());
				setEsrPipeConstraints_gen3();
				System.out.println("After esrpipeconstraints: "+problem.getConstraintsCount());
				break;
			case 5:
				addVariables_gen3();
				System.out.println("After addvariables: "+problem.getConstraintsCount());
				setObjectiveCost_esr_gen3();
				System.out.println("After setobjective: "+problem.getConstraintsCount());
				setPipeConstraints_gen3();
				System.out.println("After pipeconstraints: "+problem.getConstraintsCount());
				setHeadLossConstraints_esr_gen3();
				System.out.println("After headlossconstraints: "+problem.getConstraintsCount());
				setEsrConstraints_gen4();
				System.out.println("After esrconstraints: "+problem.getConstraintsCount());
				setEsrPipeConstraints_gen3();
				System.out.println("After esrpipeconstraints: "+problem.getConstraintsCount());
				break;
			case 6:
				addVariables_gen3();
				System.out.println("After addvariables: "+problem.getConstraintsCount());
				setObjectiveCost_esr_gen3();
				System.out.println("After setobjective: "+problem.getConstraintsCount());
				setPipeConstraints_gen3();
				System.out.println("After pipeconstraints: "+problem.getConstraintsCount());
				setHeadLossConstraints_esr_gen3();
				System.out.println("After headlossconstraints: "+problem.getConstraintsCount());
				setEsrConstraints_gen5();
				System.out.println("After esrconstraints: "+problem.getConstraintsCount());
				setEsrPipeConstraints_gen3();
				System.out.println("After esrpipeconstraints: "+problem.getConstraintsCount());
				break;
			case 7:
				addVariables_gen3();
				System.out.println("After addvariables: "+problem.getConstraintsCount());
				setObjectiveCost_esr_gen3();
				System.out.println("After setobjective: "+problem.getConstraintsCount());
				setEsrOptionConstraints();
				System.out.println("After setesroption: "+problem.getConstraintsCount());
				setPipeConstraints_gen3();
				System.out.println("After pipeconstraints: "+problem.getConstraintsCount());
				setHeadLossConstraints_esr_gen3();
				System.out.println("After headlossconstraints: "+problem.getConstraintsCount());
				setEsrConstraints_gen6();
				System.out.println("After esrconstraints: "+problem.getConstraintsCount());
				setEsrPipeConstraints_gen3();
				System.out.println("After esrpipeconstraints: "+problem.getConstraintsCount());
				break;
			case 8:
				addVariables_gen4();
				System.out.println("After addvariables: "+problem.getConstraintsCount());
				setObjectiveCost_esr_gen4();
				System.out.println("After setobjective: "+problem.getConstraintsCount());
				setEsrOptionConstraints();
				System.out.println("After setesroption: "+problem.getConstraintsCount());
				setPumpOptionConstraints();
				System.out.println("After setpumpoption: "+problem.getConstraintsCount());
				setValveOptionConstraints();
				System.out.println("After setvalveoption: "+problem.getConstraintsCount());
				setPipeConstraints_gen3();
				System.out.println("After pipeconstraints: "+problem.getConstraintsCount());
				setHeadLossConstraints_esr_gen4();
				System.out.println("After headlossconstraints: "+problem.getConstraintsCount());
				setEsrConstraints_gen4();
				System.out.println("After esrconstraints: "+problem.getConstraintsCount());
				setEsrPipeConstraints_gen3();
				System.out.println("After esrpipeconstraints: "+problem.getConstraintsCount());
				setPumpConstraints();
				System.out.println("After pumpconstraints: "+problem.getConstraintsCount());
				
				break;
			case 9:
				addVariables_gen4();
				System.out.println("After addvariables: "+problem.getConstraintsCount());
				setObjectiveCost_esr_gen4();
				System.out.println("After setobjective: "+problem.getConstraintsCount());
				setEsrOptionConstraints();
				System.out.println("After setesroption: "+problem.getConstraintsCount());
				setPumpOptionConstraints();
				System.out.println("After setpumpoption: "+problem.getConstraintsCount());
				setValveOptionConstraints();
				System.out.println("After setvalveoption: "+problem.getConstraintsCount());
				setPipeConstraints_gen3();
				System.out.println("After pipeconstraints: "+problem.getConstraintsCount());
				setHeadLossConstraints_esr_gen4();
				System.out.println("After headlossconstraints: "+problem.getConstraintsCount());
				setEsrConstraints_gen9();
				System.out.println("After esrconstraints: "+problem.getConstraintsCount());
				setEsrPipeConstraints_gen3();
				System.out.println("After esrpipeconstraints: "+problem.getConstraintsCount());
				setPumpConstraints();
				System.out.println("After pumpconstraints: "+problem.getConstraintsCount());
				break;
				
			case 10:
				addVariables_gen5();
				System.out.println("After addvariables: "+problem.getConstraintsCount());
				setObjectiveCost_esr_gen5();
				System.out.println("After setobjective: "+problem.getConstraintsCount());
				setEsrOptionConstraints_gen2();
				System.out.println("After setesroption: "+problem.getConstraintsCount());
				setPumpOptionConstraints();
				System.out.println("After setpumpoption: "+problem.getConstraintsCount());
				setValveOptionConstraints();
				System.out.println("After setvalveoption: "+problem.getConstraintsCount());
				setPipeConstraints_gen3();
				System.out.println("After pipeconstraints: "+problem.getConstraintsCount());
				setHeadLossConstraints_esr_gen5();
				System.out.println("After headlossconstraints: "+problem.getConstraintsCount());
				setEsrConstraints_gen10();
				System.out.println("After esrconstraints: "+problem.getConstraintsCount());
				setEsrPipeConstraints_gen4();
				System.out.println("After esrpipeconstraints: "+problem.getConstraintsCount());
				setPumpConstraints();
				System.out.println("After pumpconstraints: "+problem.getConstraintsCount());
				
				break;
		}
				
		System.out.println("Constraint count:" + problem.getConstraintsCount());
//				for (Constraint c : problem.getConstraints()){
//					System.out.println(c);
//				}
				
				
				
		//writeLpFile("Mokhada_input");
		//System.exit(0);
		
/////////Solver solver = factory.get(); // you should use this solver only once for one problem
		
//				SolverGLPK s = (SolverGLPK) solver;
//				Hook h = new Hook() {
//					
//					@Override
//					public void call(glp_prob arg0, glp_smcp arg1, glp_iocp arg2,
//							Map<Object, Integer> arg3) {
//						arg1.setIt_lim(1);
//					}
//				};
//				s.addHook(h);
		
		ResultStatus resultStatus = problem.solve();
		//Result result = solver.solve(problem);
		Result result = new Result(problem);
		//Result result = s.solve(problem);
		//if(result==null){
		if(resultStatus==ResultStatus.NOT_SOLVED){
			long endTime = System.currentTimeMillis();
			System.out.println("Time taken: " + (endTime - startTime) + " milliseconds");
			System.out.println("Time taken: " + problem.getTimeTaken() + " milliseconds");
			System.out.println("Nodes: "+problem.getNodes());
			throw new Exception("Optimization took too long. Please reduce input network size. If using ESR option consider removing zero demand nodes.<br> To run for extended period, please download the local version from https://www.cse.iitb.ac.in/~nikhilh/jaltantra");
		}
		else if(resultStatus==ResultStatus.OPTIMAL){
			System.out.println(result.getObjectiveValue());
			
			//set water head and pipe diameters from the result
			switch(modelNumber){
				case 0:
					setHeadsFromResult(result);
					setDiametersFromResult(result);
					//System.out.println(result);
					//printResult(result);
					break;
				case 1:
				case 2:
				case 3:
					setHeadsFromResult_esr(result);
					setDiametersFromResult(result);
					//System.out.println(result);
					//printResult(result);
					break;
				case 4:
				case 5:
				case 6:
				case 7:
					setHeadsFromResult_esr(result);
					setDiametersFromResult_gen3(result);
					//printResult_gen3(result);
					break;
				case 8:
				case 9:
					setHeadsFromResult_esr(result);
					setDiametersFromResult_gen4(result);
					//printVariables(result);
					break;
				case 10:
					setHeadsFromResult_esr_gen2(result);
					setDiametersFromResult_gen4(result);
					//printVariables_gen2(result);
					break;
				
			}
			//check();
			long endTime = System.currentTimeMillis();
			System.out.println("Time taken: " + (endTime - startTime) + " milliseconds");
			System.out.println("Time taken: " + problem.getTimeTaken() + " milliseconds");
			System.out.println("Nodes: "+problem.getNodes());
			
			//generateCoordinates();
			return true;
		}
		else{
			System.out.println("Network could not be solved.");
			return false;
			}
	}

	private void setPipeVals_ParallelLink(double minVariables[], HashMap<Pipe,Integer> pipeArray, HashMap<PipeCost,Integer> commercialPipeArray) throws Exception
	{
		//set the values of the chosen pipes
	    for(Pipe pipe:pipes.values())
	    {
	    	int pipeNo=1;
	    	//pipe's flow
	    	if(minVariables[pipeArray.size()*commercialPipeArray.size() + pipeArray.get(pipe)]>0)
	    		pipe.setFlow(minVariables[pipeArray.size()*commercialPipeArray.size() + pipeArray.get(pipe)]);
	    	else
	    		pipe.setFlow(-minVariables[pipeArray.size()*commercialPipeArray.size() + pipeArray.size() + pipeArray.get(pipe)]);
	    	for(PipeCost commPipe:pipeCost)
	    	{
	    		if(Math.abs(minVariables[(pipeArray.get(pipe)*(pipeCost.size()))+commercialPipeArray.get(commPipe)])>1E-2)
	    		{
	    			if(pipeNo==1)
	    			{
	    				pipe.setChosenPipeCost(commPipe);
	    				pipe.setDiameter(commPipe.getDiameter());
	    				pipe.setRoughness(commPipe.getRoughness());
	    				//System.out.println("Added the commercial pipe of diameter: "+commPipe.getDiameter()+" to the pipe of ID: "+pipe.getPipeID()+" for length of "+variables[(pipeArray.get(pipe)*(pipeCost.size()))+commercialPipeArray.get(commPipe)]);
	    				pipeNo++;
	    			}
	    			else if(pipeNo==2)
	    			{
	    				pipe.setChosenPipeCost2(commPipe);
	    				pipe.setDiameter2(commPipe.getDiameter());
	    				pipe.setLength2(minVariables[(pipeArray.get(pipe)*(pipeCost.size()))+commercialPipeArray.get(commPipe)]);
	    				pipe.setRoughness2(commPipe.getRoughness());
	    				pipeNo++;
	    			}
	    			else
	    				throw new Exception("more than two non-zero lengths for pipe ID: "+pipe.getPipeID());
	    		}
	    	}
	    }
	    
	    System.out.println("Variable values: First come commercial pipes with respective diamter values then,pipes with their index number, followed by commercial pipe segment's lengths in order they appeared above. Finally flows in respective pipes are printed");
	    
	    System.out.println("Commercial pipe diameters in the order they will appear after that");
	    int counter;
	    for(counter=0;counter<pipeCost.size();counter++)
	    {
	    	System.out.print(findCommPipe(counter,commercialPipeArray).getDiameter()+"\t");
	    }
	    
	    System.out.println("\nLengths of commercial pipes' segments corresponding to each pipe\n");
	    for(counter=0;counter<pipes.size();counter++)
	    {
	    	System.out.println("Pipe with ID: "+findPipe(counter,pipeArray).getPipeID());
	    	for(int temp=0;temp<pipeCost.size();temp++)
	    		System.out.print(minVariables[counter*pipeCost.size()+temp]+"\t");
	    	if(minVariables[pipeArray.size()*commercialPipeArray.size() + counter]>0)
	    		System.out.println("\n\tFlow: "+minVariables[counter+pipes.size()*pipeCost.size()]+"\n");
	    	else
	    		System.out.println("\n\tFlow: "+(-minVariables[counter+pipes.size()+pipes.size()*pipeCost.size()])+"\n");
	    }
	}
	private void setNodeVals(HashMap<Node,Integer> nodeArray,double minConstraints[], HashMap<Pipe,Integer> pipeArray)
	{
		//set the value of head at each node
	    for(int i=0;i<nodeArray.size();i++)
	    {
	    	for(Node node:nodeArray.keySet())
	    		if(nodeArray.get(node)==i)
	    		{
	    			node.setHead(source.getHead()-minConstraints[i+pipeArray.size()+1]);
	    			break;
	    		}
	    }
	}

	// Make a HashMap from pipes to numbers to allow indexing.
	private HashMap<Pipe,Integer> pipesRev()
	{
		HashMap<Pipe,Integer> pipeArray= new HashMap<Pipe,Integer>();
		int i=0;
		for(Pipe pipe:pipes.values())
	    {
	    	pipeArray.put(pipe, i);
	    	i++;
	    }
		
		return pipeArray;
		
	}
	
	// Make a HashMap from nodes to numbers to allow indexing.
	private HashMap<Node,Integer> nodesRev()
	{
		HashMap<Node,Integer> nodeArray= new HashMap<Node,Integer>();
		int i=0;
	    for(Node node:nodes.values())
	    {
	    	nodeArray.put(node, i);
	    	i++;
	    }
	    return nodeArray;
	}
	
	// Make a HashMap from pipeCost to numbers to allow indexing.
	private HashMap<PipeCost,Integer> commPipes()
	{
		HashMap<PipeCost,Integer> commercialPipeArray= new HashMap<PipeCost,Integer>();
	    int i=0;//indexing var
	    for(PipeCost commPipe:pipeCost)
	    {
	    	commercialPipeArray.put(commPipe, i);
	    	i++;
	    }
	    return commercialPipeArray;
	}
	
	// Calculate supply at source node as sum of all demands in graph.
	private double calcSourceSupply()
	{
		double sourceSupply=0;
		for(Node node:nodes.values())
    	{
    		sourceSupply+=node.getDemand();
    	}
		return sourceSupply;
	}
	
	// Make a graph object that stores the pipes each node is connected to.
	private HashMap<Node,ArrayList<Pipe>> makeGraph()
	{
		HashMap<Node,ArrayList<Pipe>> graph= new HashMap<Node,ArrayList<Pipe>>();
	    for(Node node:nodes.values())
	    	graph.put(node,new ArrayList<Pipe>());
	    for(Pipe pipe: pipes.values())
	    {
	    	Node startNode=pipe.getStartNode(),endNode=pipe.getEndNode();
	    	graph.get(startNode).add(pipe);
	    	graph.get(endNode).add(pipe);	
	    }
	    
	    return graph;
	}
	
	// Build a spanning tree using prim's algorithm.
	private HashMap<Node,ArrayList<Node>> makeSpanningTree(HashMap<Node,ArrayList<Pipe>> graph)
	{
		// Currently using prim's algorithm, can use any other algorithm as well.
		
		// spanTree is the adjacency list representation of the spanning tree.
	    HashMap<Node,ArrayList<Node>> spanTree= new HashMap<Node,ArrayList<Node>>();
	    for(Node node:nodes.values())
	    	spanTree.put(node, new ArrayList<Node>());
	    
	    // Making the min spanning tree using Prim's algorithm.
	    Set<Node> mstSet = new HashSet<Node>();
	    
	    //vertexNearNode gives the node in mstSet whose distance from the current node is stored in vertexWeight.
	    HashMap<Node,Node> vertexNearNodes = new HashMap<Node,Node>();
	    
	    // vertexWeight is the least distance of nodes from the made min spanning tree, values for vertices in the mst can be anything.
	    HashMap<Node,Double> vertexWeight = new HashMap<Node,Double>();
	    
	    // Initialise every node to be at infinite distance from source vertex.
	    for(Node node: nodes.values())
	    {
	    	vertexWeight.put(node,Double.POSITIVE_INFINITY);
	    	vertexNearNodes.put(node,source);
	    }
	    vertexWeight.put(source,0.0);
	    while(mstSet.size()!=nodes.size())
	    {
	    	Node nearestNode=source;
	    	double leastVal=Double.POSITIVE_INFINITY;
	    	for(Node next: vertexWeight.keySet())
	    	{
	    		if(!mstSet.contains(next) && vertexWeight.get(next)<leastVal)
	    		{
	    			nearestNode=next;
	    			leastVal=vertexWeight.get(next);
	    			
	    		}
	    			
	    	}
	    	// nearestNode is source in first iteration, no links are added, spanTree needn't be updated.
	    	if(nearestNode!=source)
	    	{
	    		spanTree.get(nearestNode).add(vertexNearNodes.get(nearestNode));
	    		spanTree.get(vertexNearNodes.get(nearestNode)).add(nearestNode);
	    	}
	    	// vertexWeight and vertexNearNodes updated according to node to be added.
	    	for(Pipe addedPipe : graph.get(nearestNode))
			{
				//if the addedPipe starts at nearestNode
	    		if(nearestNode.getNodeID()==addedPipe.getStartNode().getNodeID())
	    		{
	    			//check if this pipe is secondary link for sure, then update the fcm matrix
	    			/*if(mstSet.contains(addedPipe.getEndNode()))
	    			{
	    				
	    			}*/
	    			
	    			if(vertexWeight.get(addedPipe.getEndNode())>addedPipe.getLength())
    				{
	    				vertexWeight.put(addedPipe.getEndNode(), addedPipe.getLength());
	    				vertexNearNodes.put(addedPipe.getEndNode(),nearestNode);
    				}
	    		}
	    		else//if addedPipe ends at nearestNode
	    		{
	    			if(vertexWeight.get(addedPipe.getStartNode())>addedPipe.getLength())
    				{
	    				vertexWeight.put(addedPipe.getStartNode(), addedPipe.getLength());
	    				vertexNearNodes.put(addedPipe.getStartNode(),nearestNode);
    				}
	    		}
			}
	    	
	    	mstSet.add(nearestNode);
	    	
	    }
	    return spanTree;
	}
	
	// Perform bfs on graph given in adjacency list form.
	// Return a 2-element list of (HashMap from node to distance from root, HashMap from node to it's parent).
	private List<Object> bfs(HashMap<Node,ArrayList<Node>> graph,Node graphSource)
	{
	    HashMap<Node,Double> dist = new HashMap<Node,Double>();
	    HashMap<Node,Node> par = new HashMap<Node,Node>();
	    for(Node node: graph.keySet())
	    {
	    	dist.put(node,-1.0);
	    	//par.put(node,source);
	    }
	    dist.put(graphSource,0.0);
	    Queue<Node> bfsQueue=new LinkedList<Node>();
	    bfsQueue.add(graphSource);
	    while(!bfsQueue.isEmpty())
	    {
	    	Node temp=bfsQueue.remove();
	    	for(Node node:graph.get(temp))
	    	{
	    		if(dist.get(node)==-1)
	    		{
	    			dist.put(node,dist.get(temp)+1);
	    			par.put(node,temp);
	    			bfsQueue.add(node);
	    		}
	    	}
	    }
	    return Arrays.asList(dist,par);
	}
	
	// Make matrix C in the paper to tell direction of consideration of each pipe in each cycle 
	// for headloss in loop constraint calculation.
	private int[][] makeCFDMat(HashMap<Node,ArrayList<Pipe>> graph,HashMap<Node,ArrayList<Node>> spanTree,HashMap<Node,Double> dist,HashMap<Node,Node> par,HashMap<Pipe,Integer> pipeArray)
	{
		int cycleFlowDir[][]= new int [pipes.size()-nodes.size()+1][pipes.size()];
	    int cycleNo=0;
	    for(Pipe pipe: pipes.values())
	    {
	    	Node startNode=pipe.getStartNode(),endNode=pipe.getEndNode();
	    	Node node1=startNode,node2=endNode;
	    	if(!spanTree.get(startNode).contains(endNode))//secondary link
    		{
	    		cycleFlowDir[cycleNo][pipeArray.get(pipe)]=1;
	    		
	    		if(dist.get(startNode)>dist.get(endNode))
	    		{
	    			Node parent=par.get(startNode);
	    			while(dist.get(node1)>dist.get(endNode))
	    			{
	    				for(Pipe possPipe:graph.get(node1))
	    				{
	    					if(possPipe.getStartNode().getNodeID()==parent.getNodeID())
							{
	    						cycleFlowDir[cycleNo][pipeArray.get(possPipe)]=1;
	    						break;
							}
	    					else if(possPipe.getEndNode().getNodeID()==parent.getNodeID())
	    					{
	    						cycleFlowDir[cycleNo][pipeArray.get(possPipe)]=-1;
	    						break;
	    					}
	    						
	    				}
	    				node1=parent;
	    				parent=par.get(parent);
	    			}
	    		}
	    		else if(dist.get(startNode)<dist.get(endNode))
	    		{
	    			Node parent=par.get(endNode);
	    			while(dist.get(node2)>dist.get(startNode))
	    			{
	    				for(Pipe possPipe:graph.get(node2))
	    				{
	    					if(possPipe.getStartNode().getNodeID()==parent.getNodeID())
	    					{
	    						cycleFlowDir[cycleNo][pipeArray.get(possPipe)]=-1;
	    						break;
	    					}
	    					else if(possPipe.getEndNode().getNodeID()==parent.getNodeID())
	    					{
	    						cycleFlowDir[cycleNo][pipeArray.get(possPipe)]=1;
	    						break;
	    					}
	    				}
	    				node2=parent;
	    				parent=par.get(parent);
	    			}
	    		}
	    		while(node1.getNodeID()!=node2.getNodeID())
	    		{
	    			Node parent1=par.get(node1),parent2=par.get(node2);
	    			for(Pipe possPipe:graph.get(node1))
					{
						if(possPipe.getStartNode().getNodeID()==parent1.getNodeID())
						{
							cycleFlowDir[cycleNo][pipeArray.get(possPipe)]=1;
							break;
						}
						else if(possPipe.getEndNode().getNodeID()==parent1.getNodeID())
						{
							cycleFlowDir[cycleNo][pipeArray.get(possPipe)]=-1;
							break;
						}
							
					}
	    			node1=parent1;
	    			parent1=par.get(parent1);
	    			for(Pipe possPipe:graph.get(node2))
					{
						if(possPipe.getStartNode().getNodeID()==parent2.getNodeID())
						{
							cycleFlowDir[cycleNo][pipeArray.get(possPipe)]=-1;
							break;
						}
						else if(possPipe.getEndNode().getNodeID()==parent2.getNodeID())
						{
							cycleFlowDir[cycleNo][pipeArray.get(possPipe)]=1;
							break;
						}
					}
					node2=parent2;
					parent2=par.get(parent2);
	    		}
	    		cycleNo++;
    		}
	    }
	    return cycleFlowDir;
	}

	// Make the matrix S to tell direction of consideration of pipes while reaching each node
    // for headloss to reach a node calculations.
	private int[][] makePLFMat(HashMap<Node,Node> par,HashMap<Node,ArrayList<Pipe>> graph,HashMap<Pipe,Integer> pipeArray,HashMap<Node,Integer> nodeArray)
	{
		int primLinkFlow[][]= new int[nodes.size()][pipes.size()];
	    for(Node node:nodes.values())
	    {
	    	Node temp=node;
	    	while(temp.getNodeID()!=source.getNodeID())
	    	{
	    		Node parent=par.get(temp);
	    		for(Pipe possPipe:graph.get(temp))
				{
					if(possPipe.getStartNode().getNodeID()==parent.getNodeID())
					{
						primLinkFlow[nodeArray.get(node)][pipeArray.get(possPipe)]=1;
						break;
					}
					else if(possPipe.getEndNode().getNodeID()==parent.getNodeID())
					{
						primLinkFlow[nodeArray.get(node)][pipeArray.get(possPipe)]=-1;
						break;
					}
						
				}
	    		temp=parent;
	    	}
	    }
	    return primLinkFlow;
	}
	
	// Make the matrix F to tell whether a pipe starts or ends at a node
    // for total flow to a node calculations.
	private int[][] makeFDMat(HashMap<Pipe,Integer> pipeArray,HashMap<Node,Integer> nodeArray)
	{
		int flowDir[][]= new int[nodes.size()][pipes.size()];
	    for(Pipe pipe:pipes.values())
	    {
	    	flowDir[nodeArray.get(pipe.getStartNode())][pipeArray.get(pipe)]=-1;
	    	flowDir[nodeArray.get(pipe.getEndNode())][pipeArray.get(pipe)]=1;
	    }
	    return flowDir;
	}
	private void setOptions(Cycle_Optimizer cycleOptimizer)
	{
		//System.out.println(cycleOptimizer.setStringOption("print_info_string","yes"));
	    //System.out.println(cycleOptimizer.setStringOption("warm_start_init_point","yes"));
	    //System.out.println(cycleOptimizer.setStringOption("nlp_scaling_method","none"));
	    System.out.println(cycleOptimizer.setIntegerOption("max_iter",1000000));
	    System.out.println(cycleOptimizer.setNumericOption("constr_viol_tol",0.00001));
	    //System.out.println(cycleOptimizer.setNumericOption("acceptable_tol",0.001));
//	    System.out.println(cycleOptimizer.setStringOption("derivative_test","second-order"));
	    //System.out.println(cycleOptimizer.setStringOption("derivative_test_print_all","yes"));
	    //System.out.println(cycleOptimizer.setNumericOption("max_cpu_time",60));
	    //System.out.println(cycleOptimizer.setNumericOption("point_perturbation_radius",0));
	    //System.out.println(cycleOptimizer.setIntegerOption("print_level",12));
	    //System.out.println(cycleOptimizer.setStringOption("output_file","opfile.txt"));
	    //System.out.println(cycleOptimizer.setIntegerOption("file_print_level",6));
	}
	
	private void setPipeVals_DiscSegment(double variables[], HashMap<Pipe,Integer> pipeArray, HashMap<PipeCost,Integer> commercialPipeArray) throws Exception
	{
		//set the values of the chosen pipes
	    for(Pipe pipe:pipes.values())
	    {
	    	int pipeNo=1;
	    	//pipe's flow
	    	pipe.setFlow(variables[pipeArray.size()*commercialPipeArray.size() + pipeArray.get(pipe)]);
	    	
	    	for(PipeCost commPipe:pipeCost)
	    	{
	    		if(Math.abs(variables[(pipeArray.get(pipe)*(pipeCost.size()))+commercialPipeArray.get(commPipe)])>1E-2)
	    		{
	    			if(pipeNo==1)
	    			{
	    				pipe.setChosenPipeCost(commPipe);
	    				pipe.setDiameter(commPipe.getDiameter());
	    				pipe.setRoughness(commPipe.getRoughness());
	    				//System.out.println("Added the commercial pipe of diameter: "+commPipe.getDiameter()+" to the pipe of ID: "+pipe.getPipeID()+" for length of "+variables[(pipeArray.get(pipe)*(pipeCost.size()))+commercialPipeArray.get(commPipe)]);
	    				pipeNo++;
	    			}
	    			else if(pipeNo==2)
	    			{
	    				pipe.setChosenPipeCost2(commPipe);
	    				pipe.setDiameter2(commPipe.getDiameter());
	    				pipe.setLength2(variables[(pipeArray.get(pipe)*(pipeCost.size()))+commercialPipeArray.get(commPipe)]);
	    				pipe.setRoughness2(commPipe.getRoughness());
	    				pipeNo++;
	    			}
	    			else
	    				throw new Exception("more than two non-zero lengths for pipe ID: "+pipe.getPipeID());
	    		}
	    	}
	    }
	    System.out.println("Variable values: First come commercial pipes with respective diamter values then,pipes with their index number, followed by commercial pipe segment's lengths in order they appeared above. Finally flows in respective pipes are printed");
	    
	    System.out.println("Commercial pipe diameters in the order they will appear after that");
	    int counter;
	    for(counter=0;counter<pipeCost.size();counter++)
	    {
	    	System.out.print(findCommPipe(counter,commercialPipeArray).getDiameter()+"\t");
	    }
	    
	    System.out.println("\nLengths of commercial pipes' segments corresponding to each pipe\n");
	    for(counter=0;counter<pipes.size();counter++)
	    {
	    	System.out.println("Pipe with ID: "+findPipe(counter,pipeArray).getPipeID());
	    	for(int temp=0;temp<pipeCost.size();temp++)
	    		System.out.print(variables[counter*pipeCost.size()+temp]+"\t");
	    	System.out.println("\n\tFlow: "+variables[counter+pipes.size()*pipeCost.size()]+"\n");
	    }
	}
	
	/**return all orientations of given graph where each orientation
	 * maps each pipe to +1,-1 so that -1 or +1 means towards or away from start of pipe*/
	private ArrayList<HashMap<Pipe,Integer>> getOrientations(HashMap<Node,ArrayList<Pipe>> graph,HashMap<Node,Integer> nodeArray,HashMap<Pipe,Integer> pipeArray,int maxOrientAtEachStep)
	{
//		HashMap<Integer,ArrayList<Integer>> adjList=new HashMap<Integer,ArrayList<Integer>>();
//		HashMap<Node,Boolean> visit=new HashMap<Node,Boolean>();
////		HashMap<Integer,Set<Integer>> edgeMapping=new HashMap<Integer,Set<Integer>>();
//		HashMap<Node,Integer> level=new HashMap<Node,Integer>();
//		Set<Pipe> remEdge=new HashSet<Pipe>();
//		
//		for(Node n:nodes.values())
//		{
//			ArrayList<Integer> newarli=new ArrayList<Integer>();
//			adjList.put(nodeArray.get(n),newarli);
//			visit.put(n,false);
//			level.put(n,-1);
//		}
//		level.put(source, 0);
//		
//		int bla=redGraphDFS(graph,nodeArray,pipeArray,adjList,0,source,visit,edgeMapping,false,0,level,remEdge,bridges,pipeStart,pipeEnd);
//		return;
		ArrayList<Integer> sourceList=new ArrayList<Integer>();
		sourceList.add(nodeArray.get(source));
		redGraph rootGraph=new redGraph(nodeArray.get(source),sourceList,maxOrientAtEachStep);
		for(Node n:nodes.values())
		{
			ArrayList<edgeGroup> newarli=new ArrayList<edgeGroup>();
			rootGraph.adjList.put(nodeArray.get(n),newarli);
		}
		for(Pipe p:pipes.values())
		{
			int node2=nodeArray.get(p.getEndNode()),node1=nodeArray.get(p.getStartNode());
			edgeGroup pg=new edgeGroup(pipeArray.get(p),node1,node2);
			rootGraph.adjList.get(node1).add(pg);
			rootGraph.adjList.get(node2).add(pg);
		}
		
		redGraph endGraph=rootGraph.reduce();
		ArrayList<HashMap<edgeGroup,Integer>> orientations=endGraph.enumerateOrientations();
		return egToPipe(orientations,pipeArray);
		
		
	}
	
	
	private ArrayList<HashMap<Pipe,Integer>> getOrientations1(HashMap<Node,ArrayList<Pipe>> graph,HashMap<Node,Integer> nodeArray,HashMap<Pipe,Integer> pipeArray, int maxNoOrient)
	{
		ArrayList<HashMap<Pipe,Integer>> orientations=new ArrayList<HashMap<Pipe,Integer>>();
		File file=new File("/Users/saumyagoyal/Desktop/or.txt");
		try
		{
			Scanner sc = new Scanner(file);
			int lines=0;
			while(sc.hasNextInt())
			{
				lines+=1;
				System.out.println("Reading line no: "+lines);
				HashMap<Pipe,Integer> curr_orient=new HashMap<Pipe,Integer>();
				for(int i=0;i<pipeArray.size();i++)
				{
					int id=sc.nextInt();
					int or=sc.nextInt();
					if(or==0)
						System.out.println("Found 0 for id: "+id);
					curr_orient.put(pipes.get(id), or);
				}
				System.out.println("Line read");
				orientations.add(curr_orient);
			}
			System.out.println("File read, returning orientations to be used");
			sc.close();
			return orientations;
		}
		catch(Exception e)
		{
			System.out.println("Exception: "+e);
			System.out.println("Calling the earlier version of getOrientations with 10 orientations");
			return getOrientations(graph,nodeArray,pipeArray,maxNoOrient);
		}
	}
	
	/**convert orientations from edgeGroup to pipe object*/
	private ArrayList<HashMap<Pipe,Integer>> egToPipe(ArrayList<HashMap<edgeGroup,Integer>> egOrientations,HashMap<Pipe,Integer> pipeArray)
	{
		ArrayList<HashMap<Pipe,Integer>> orientations=new ArrayList<HashMap<Pipe,Integer>>();
		for(HashMap<edgeGroup,Integer> orient:egOrientations)
		{
			HashMap<Pipe,Integer> pOrientation=new HashMap<Pipe,Integer>();
			for(edgeGroup p:orient.keySet())
			{
				//orient.get(p) works as start node of pipe is node1 of edgegroup and end node is node2
				pOrientation.put(findPipe(p.getID(),pipeArray), orient.get(p));
			}
			orientations.add(pOrientation);
		}
		return orientations;
		
	}
	
	private boolean cyclicOptimize(long startTime, Formulation formulation) throws Exception
	{
		
		final HashMap<Pipe, Integer> pipeArray = pipesRev();
	    final HashMap<Node, Integer> nodeArray = nodesRev();
	    final HashMap<PipeCost, Integer> commercialPipeArray = commPipes();
	    
	    final double minFlow = generalProperties.min_flow;
	    final double sourceSupply = calcSourceSupply();

	    final double maxFlow;
	    if(generalProperties.max_flow == 0) {
	    	// Maximum flow is unspecified.
	    	maxFlow = sourceSupply;
	    }
	    else
	    	maxFlow = Math.min(generalProperties.max_flow, sourceSupply);
	    source.setDemand(-sourceSupply);
	    
	    
	    // Graph stores the pipes each node is connected to.
	    final HashMap<Node,ArrayList<Pipe>> graph= makeGraph();
//	    final Set<Pipe> bridges=bridgeSet(graph);
//	    HashMap<Integer,ArrayList<Integer>> redAdjList=new HashMap<Integer,ArrayList<Integer>>();
//	    HashMap<Integer,Set<Integer>> redEdgeMapping=new HashMap<Integer,Set<Integer>>();
//	    HashMap<Integer,Integer> pipeStart=new HashMap<Integer,Integer>(),pipeEnd=new HashMap<Integer,Integer>();
//	    reduceGraph(graph,nodeArray,pipeArray,bridges,redAdjList,redEdgeMapping,pipeStart,pipeEnd);
	    
	    // spanTree is adjacency list representation of the spanning tree.
	    final HashMap<Node,ArrayList<Node>> spanTree= makeSpanningTree(graph);
	    
	    
	    // Carrying out BFS on obtained spanning tree.
	    final List<Object> bfsRes=bfs(spanTree,source);
	    final HashMap<Node,Double> dist = (HashMap<Node,Double>)bfsRes.get(0);
	    final HashMap<Node,Node> par = (HashMap<Node,Node>)bfsRes.get(1);
	    
	    // Matrix C in paper.
	    final int cycleFlowDir[][]= makeCFDMat(graph,spanTree,dist,par,pipeArray);
	    
	    // Matrix S in paper.
	    final int primLinkFlow[][]= makePLFMat(par,graph,pipeArray,nodeArray);
	    
	    // Matrix F in paper.
	    final int flowDir[][]= makeFDMat(pipeArray,nodeArray);
	    
	    final int no_pipes=pipes.size(),no_nodes=nodes.size(),no_cycles,no_commPipes=pipeCost.size();
	    no_cycles = no_pipes-no_nodes+1;
	    
	    Cycle_Optimizer cycleOptimizer = new Cycle_Optimizer(pipeArray,nodeArray,commercialPipeArray,flowDir,cycleFlowDir,primLinkFlow,minFlow,maxFlow,source.getHead(),pipes, formulation);

	    setOptions(cycleOptimizer);

	    int status;
	    Boolean solved=false;
	    double minCost=Double.MAX_VALUE,minVariables[]=new double[0],minConstraints[]=new double[0];

	    int maxNoOrient = 10;
	    long solvingStartTime=System.currentTimeMillis();

	    ArrayList<HashMap<Pipe,Integer>> usedOrientations=new ArrayList<HashMap<Pipe,Integer>>();
	    ArrayList<HashMap<Pipe,Integer>> finalOrientations=new ArrayList<HashMap<Pipe,Integer>>();

	    if (formulation == Formulation.RE_SOLVING){
	    	//each element of list is an orientation
		    //each orientation maps each pipe to +1,-1 so that -1 or +1 means towards or away from start of pipe
		    ArrayList<HashMap<Pipe,Integer>> orientations = getOrientations1(graph,nodeArray,pipeArray,maxNoOrient);

		    int noOp = 1,incr = 1;
		    int noOptPerOrientation=1;
		    solvingStartTime=System.currentTimeMillis();
		    System.out.println("Number of orientations: "+orientations.size());
		    System.out.println("Will do "+noOptPerOrientation+" optimisations per orientation");
		    incr=Math.max(orientations.size()/maxNoOrient,1);
		    for(int orientationNo=0;orientationNo<orientations.size();orientationNo=orientationNo+incr)
		    {
		    	HashMap<Pipe,Integer> or=orientations.get(orientationNo);
		    	for(int or_opt_num=0;or_opt_num<noOptPerOrientation;or_opt_num++)
		    	{
			    	System.out.println("Doing "+noOp+"th optimisation");
			    	long sTime=System.currentTimeMillis();
			    	status=cycleOptimizer.OptimizeNLP(or);
			    	
			    	//proper condition for solve successful/unsuccessful not added yet
			    	System.out.println("Done "+noOp+"th optimisation");
			    	noOp++;
				    double cost=cycleOptimizer.getObjectiveValue(),variables[]=cycleOptimizer.getVariableValues(),constraints[]=cycleOptimizer.getConstraintValues();
				    System.out.println("Optimized cost is: "+cost);
				    long endTime = System.currentTimeMillis();
				    System.out.println("Time taken for current optimisation is: "+(endTime-sTime)+" milliseconds");
				    usedOrientations.add(or);
				    HashMap<Pipe,Integer> thisFinalOr=new HashMap<Pipe,Integer>();
				    for(Pipe p:pipeArray.keySet())
				    {
				    	thisFinalOr.put(p, (int)Math.signum(variables[pipeArray.get(p)+no_pipes*no_commPipes]));
				    }
				    finalOrientations.add(thisFinalOr);
				    if(status==cycleOptimizer.SOLVE_SUCCEEDED)
				    {
			    		solved=true;
					    if(cost<minCost)
					    {
					    	minVariables=variables.clone();
					    	minConstraints=constraints.clone();
						    minCost=cost;
					    }
				    }
				    int minutesGone=(int)((endTime-solvingStartTime)/(1000*60));
				    System.out.println("Currently going at "+(minutesGone/(noOp-1))+"min per orientation");
		    	}
		    }
		    printOrientations(usedOrientations,finalOrientations, formulation);
		    //set the values of the chosen pipes
		    setPipeVals_DiscSegment(minVariables,pipeArray,commercialPipeArray);
		    //set the value of head at each node
		    setNodeVals(nodeArray,minConstraints,pipeArray);
		    System.out.println("\n\n\n");
		    System.out.println("Least cost is: "+minCost);
		    long endTime = System.currentTimeMillis();
		    System.out.println("Total time taken for "+(noOp-1)+" optimisations is: "+(endTime-startTime)+" millisecond");
	    }
	    else {
	    	int noOp = maxNoOrient;
	    	System.out.println("Will be doing " + maxNoOrient + " orientations");
	    	for(;noOp>0;noOp--)
		    {
		    	HashMap<Pipe,Integer> finalOrient=new HashMap<Pipe,Integer>();
		    	System.out.println("Doing "+(maxNoOrient + 1 - noOp)+"th optimisation");
		    	long sTime=System.currentTimeMillis();
		    	status=cycleOptimizer.OptimizeNLP(maxNoOrient + 1 - noOp);
		    	
		    	System.out.println("Done "+(maxNoOrient + 1 - noOp)+"th optimisation");
			    double cost=cycleOptimizer.getObjectiveValue(),variables[]=cycleOptimizer.getVariableValues(),constraints[]=cycleOptimizer.getConstraintValues();
			    System.out.println("Optimized cost is: "+cost);
			    long endTime = System.currentTimeMillis();
			    System.out.println("Time taken for current optimisation is: "+(endTime-sTime)+" milliseconds");
			    usedOrientations.add(new HashMap<Pipe,Integer>(cycleOptimizer.getOrientation()));
			    
			    for(Pipe p:pipeArray.keySet())
			    {
			    	if (formulation == Formulation.PARALLEL_LINK)
			    	{// PL formulation.
			    		if(variables[pipeArray.size()*commercialPipeArray.size() + pipeArray.get(p)]>1E-2)
				    		finalOrient.put(p, 1);
				    	else
				    		finalOrient.put(p, -1);
			    	}
			    	else
			    	{
			    		finalOrient.put(p, (int)Math.signum(variables[pipeArray.size()*commercialPipeArray.size() + pipeArray.get(p)]));
			    	}
			    }
			    
			    finalOrientations.add(finalOrient);
			    if(status==cycleOptimizer.SOLVE_SUCCEEDED)
			    {
		    		solved=true;
				    if(cost<minCost)
				    {
				    	minVariables=variables.clone();
				    	minConstraints=constraints.clone();
					    minCost=cost;
				    }
			    }
			    int minutesGone=(int)((endTime-solvingStartTime)/(1000*60));
				// TODO: Verify if the below statement will print the correct value or not
			    System.out.println("Currently going at "+(minutesGone/(maxNoOrient+1-noOp))+"min per orientation");
		    }

		    printOrientations(usedOrientations,finalOrientations, formulation);
	    
		    System.out.println("\n\n\n");
		    System.out.println("Least cost is: "+minCost);
		    long endTime = System.currentTimeMillis();
			// TODO: Verify if the below statement will print the correct value or not
			System.out.println("Total time taken for "+noOp+" optimisations is: "+(endTime-startTime)+" millisecond");
		    
		    //set the values of the chosen pipes
		    if (formulation == Formulation.PARALLEL_LINK) {
		    	setPipeVals_ParallelLink(minVariables,pipeArray,commercialPipeArray);
		    }
		    else {
		    	setPipeVals_DiscSegment(minVariables,pipeArray,commercialPipeArray);
		    }
		    //set the value of head at each node
		    setNodeVals(nodeArray,minConstraints,pipeArray);
	    }

	    cycleOptimizer.dispose();
	    return solved;
	}
	
	// Print the starting point and final point orientations.
	private void printOrientations(ArrayList<HashMap<Pipe,Integer>> orientations,ArrayList<HashMap<Pipe,Integer>> finalOrientations, Formulation formulation)
	{
		if (formulation == Formulation.RE_SOLVING) {
			System.out.println("Next "+orientations.size()+" lines represent all the input orientations in the format of pipeid orientation pipeid orientation ...");
		}
		else{
			System.out.println("The following "+orientations.size()+" lines represent the starting point orientations used in the order of pipeid orientation pipeid orientation ...");
		}
		for(HashMap<Pipe,Integer> or:orientations)
		{
			for(Pipe p:or.keySet())
			{
				System.out.print(p.getPipeID()+"\t"+or.get(p)+"\t");
			}
			System.out.println("");
		}

		if (formulation == Formulation.RE_SOLVING) {
			System.out.println("Next "+finalOrientations.size()+" lines represent all the final orientations in the format of pipeid orientation pipeid orientation ...");
		}
		else {
			System.out.println("The following "+finalOrientations.size()+" lines represent the final point orientations used in the order of pipeid orientation pipeid orientation ...");
		}
		for(HashMap<Pipe,Integer> or:finalOrientations)
		{
			for(Pipe p:or.keySet())
			{
				System.out.print(p.getPipeID()+"\t"+or.get(p)+"\t");
			}
			System.out.println("");
		}
	}

	private Node findNode(int index,HashMap<Node,Integer> nodeArray)
	{
		for(Node node:nodeArray.keySet())
			if(nodeArray.get(node)==index)
				return node;
		return findNode(0,nodeArray);
	}
	private Pipe findPipe(int index,HashMap<Pipe,Integer> pipeArray)
	{
		for(Pipe pipe:pipeArray.keySet())
			if(pipeArray.get(pipe)==index)
				return pipe;
		return findPipe(0,pipeArray);
	}
	private PipeCost findCommPipe(int index,HashMap<PipeCost,Integer> commercialPipeArray)
	{
		for(PipeCost commPipe:commercialPipeArray.keySet())
			if(commercialPipeArray.get(commPipe)==index)
				return commPipe;
		return findCommPipe(0,commercialPipeArray);
	}
	

	//generate coordinate string to use for EPANET output file generation
	private void generateCoordinatesHelper(Node head, double x, double y, double minx, double maxx){
		coordinatesString += head.getNodeID()+" "+x+" "+y+",";
		List<Pipe> outgoingPipes = head.getOutgoingPipes();
		int n = outgoingPipes.size();
		
		if(n>1){
			double min_length = Double.MAX_VALUE;
			for(Pipe p : outgoingPipes){
				min_length = Math.min(min_length, p.getLength());
			}
			double max_width = min_length/(1-(1/(double)n));
			minx = Math.max(minx, x-max_width);
			maxx = Math.min(maxx, x+max_width);
		}
		for(int i=0;i<n;i++){
			Pipe p = outgoingPipes.get(i);
			double newminx = minx + (maxx-minx)*i/n;
			double newmaxx = minx + (maxx-minx)*(i+1)/n;
			double newx = (newminx+newmaxx)/2;
			double newy = y + Math.sqrt(p.getLength()*p.getLength() - (newx-x)*(newx-x));
			generateCoordinatesHelper(p.getEndNode(), newx, newy, newminx, newmaxx);
		}
		
	}

	public HashMap<Integer,Node> getNodes(){
		return nodes;
	}
	
	public HashMap<Integer,Pipe> getPipes(){
		return pipes;
	}
	
	public List<PipeCost> getPipeCost(){
		return pipeCost;
	}
		
}

