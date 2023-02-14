package optimizer;

import java.util.*;
import org.coinor.Ipopt;
import optimizer.Formulation;

//n=no of var, m=no of constraints
public class Cycle_Optimizer extends Ipopt
{
	final int n,m,count_bounds=0,dcount_start=0;//rename the vars
	final int no_nodes,no_pipes,no_cycles,no_commPipes;
	final double minFlow,maxFlow,sourceHead,lambda=1.216*Math.pow(10, -4),gamma=1.209*Math.pow(10,10);
	final HashMap<Pipe,Integer> pipeArray;
	final HashMap<Node,Integer> nodeArray;
	final HashMap<PipeCost,Integer> commPipeArray;
	final int flowDir[][],cycleFlowDir[][],primLinkFlow[][];
	double start=0.5;
	long seed=-1;
	final List<Pipe> inputPipes;
	HashMap<Pipe,Integer> orientation=null;

	final Formulation formulation;
	
	public Cycle_Optimizer(HashMap<Pipe, Integer> pipeArray, HashMap<Node, Integer> nodeArray,
			HashMap<PipeCost,Integer> commPipeArray, int flowDir[][],
			int cycleFlowDir[][], int primLinkFlow[][], double minFlow,
			double maxFlow, double sourceHead, HashMap<Integer, Pipe> inputPipes, Formulation formulation)
	{
		no_nodes=nodeArray.size();
		no_pipes=pipeArray.size();
		no_commPipes=commPipeArray.size();
		
		int nele_jac, nele_hess;
		
		if (formulation == Formulation.PARALLEL_LINK) {
			this.n = (no_commPipes+2)*no_pipes;
			this.m = no_nodes+3*no_pipes+1;
			nele_jac = no_pipes*(2*no_nodes+(no_pipes+2)*(no_commPipes+2));
			nele_hess = (3+2*no_commPipes)*no_pipes;
		} else {
			this.n = (no_commPipes+1)*no_pipes;
			this.m = no_nodes+3*no_pipes+1;
			nele_jac = no_pipes*(no_pipes+no_nodes+2)*(1+no_commPipes);
			nele_hess = (1+no_commPipes)*no_pipes;
		}
		this.pipeArray=pipeArray;
		this.nodeArray=nodeArray;
		this.commPipeArray=commPipeArray;
		
		this.flowDir=flowDir;
		this.cycleFlowDir=cycleFlowDir;
		this.primLinkFlow=primLinkFlow;
		this.minFlow=minFlow;
		this.maxFlow=maxFlow;
		this.sourceHead=sourceHead;
		this.no_cycles=this.no_pipes-this.no_nodes+1;

		if (formulation != Formulation.RE_SOLVING) {
			orientation = new HashMap<Pipe,Integer>();
		}
		//deterministic ordering of pipes based on input file
		this.inputPipes=new ArrayList<Pipe>();
		List<Integer> tempIndList=new ArrayList<Integer>(inputPipes.keySet());
		Collections.sort(tempIndList);
		for(Integer i:tempIndList)
		{
			this.inputPipes.add(inputPipes.get(i));
		}
		this.formulation = formulation;
		create(n,m,nele_jac,nele_hess,Ipopt.C_STYLE);
//		System.out.println("Object created successfully");
	}
	
	public HashMap<Pipe,Integer> getOrientation()
	{
		return orientation;
	}

	// Return Node object corresponding to given index.
	private Node findNode(int index)
	{
		for(Node node:nodeArray.keySet())
			if(nodeArray.get(node)==index)
				return node;
		return findNode(0);
	}
	
	// Return Pipe object corresponding to given index.
	private Pipe findPipe(int index)
	{
		for(Pipe pipe:pipeArray.keySet())
			if(pipeArray.get(pipe)==index)
				return pipe;
		return findPipe(0);
	}
	
	// Return PipeCost object corresponding to given index.
	private PipeCost findCommPipe(int index)
	{
		for(PipeCost commPipe:commPipeArray.keySet())
			if(commPipeArray.get(commPipe)==index)
				return commPipe;
		return findCommPipe(0);
	}
	
	
	protected boolean get_bounds_info(int n, double[] x_L, double[] x_U, int m, double[] g_L, double[] g_U)
	{
//		System.out.println("Entering get_bounds_info");
		int i,j;

		// Variable bounds for length segment variables.
		for(i=0;i<no_pipes;i++)
		{
			//considering link i
			for(j=0;j<no_commPipes;j++)
			{
				//length of segment containing jth commercial pipe
				x_L[i*no_commPipes+j]=0;
				x_U[i*no_commPipes+j]=findPipe(i).getLength();
			}
		}

		// Variable bounds for flow variables.
		if (formulation == Formulation.PARALLEL_LINK) {
			for(i=0;i<no_pipes;i++)
			{
				x_L[i+no_pipes*no_commPipes]=minFlow;
				x_U[i+no_pipes*no_commPipes]=maxFlow;
				x_L[i+no_pipes+no_pipes*no_commPipes]=minFlow;
				x_U[i+no_pipes+no_pipes*no_commPipes]=maxFlow;
			}
		} else if (formulation == Formulation.DISCRETE_SEGMENT) {
			for(i=0;i<no_pipes;i++) {
				x_L[i+no_pipes*no_commPipes]=-maxFlow;
				x_U[i+no_pipes*no_commPipes]=maxFlow;
			}
		}
		else {
			assert(orientation != null);
			for(Pipe p:orientation.keySet())
			{
				if(orientation.get(p)==1)
				{
					x_L[pipeArray.get(p)+no_pipes*no_commPipes]=minFlow;
					x_U[pipeArray.get(p)+no_pipes*no_commPipes]=maxFlow;
				}
				else if(orientation.get(p)==-1)
				{
					x_L[pipeArray.get(p)+no_pipes*no_commPipes]=-maxFlow;
					x_U[pipeArray.get(p)+no_pipes*no_commPipes]=-minFlow;
				}
				else
				{
					x_L[pipeArray.get(p)+no_pipes*no_commPipes]=-maxFlow;
					x_U[pipeArray.get(p)+no_pipes*no_commPipes]=maxFlow;
				}
			}

		}
		
		// Constraint bounds for total flow through a node.
		for(i=0;i<no_nodes;i++)
		{
			g_L[i]=0;
			g_U[i]=0;
		}
		// Constraint bounds for pressure reduction in each cycle.
		for(i=0; i <no_cycles; i++ ) 
		{
			g_L[i+no_nodes] = 0;
			g_U[i+no_nodes] = 0;
		}
		// Constraint bounds for pressure reduction till node is reached.
		for(i = 0; i < no_nodes; i++ ) 
		{
			Node currNode=findNode(i);
			g_L[i+no_cycles+no_nodes] = 0;
			g_U[i+no_nodes+no_cycles] = sourceHead-currNode.getElevation()-currNode.getResidualPressure();
		}
		
		// Constraint bounds for sensible value of flow.
		if (formulation == Formulation.PARALLEL_LINK) {
			for( i = 0; i < no_pipes; i++ ) 
			{
				g_L[i+no_cycles+no_nodes*2] = 0;
				g_U[i+no_cycles+no_nodes*2] = 0;
			}
		}
		else {
			for( i = 0; i < no_pipes; i++ ) 
			{
				g_L[i+no_cycles+no_nodes*2] = minFlow;
				g_U[i+no_cycles+no_nodes*2] = maxFlow;
			}
		}
		
		// Constraint bounds for sum of lengths of pipe segments in each link.
		for( i = 0; i < no_pipes; i++ ) 
		{
			g_L[i+no_cycles+no_nodes*2+no_pipes] = findPipe(i).getLength();
			g_U[i+no_cycles+no_nodes*2+no_pipes] = findPipe(i).getLength();
		}
		return true;
	}
	
	protected boolean get_starting_point(int n, boolean init_x, double[] x, boolean init_z, double[] z_L, double[] z_U, int m, boolean init_lambda, double[] lambda) 
	{
		assert init_z == false;
		assert init_lambda = false;
		java.util.Random generator=new Random();
		
		if(this.seed!=-1)
		{// set the seed
			generator.setSeed(this.seed);
			System.out.println("Seeded the value "+this.seed);
		}
		
		//if( init_x ) 
		{
			// Starting point for lenths of each diamter.
			for(Pipe p:this.inputPipes)
			{
				int i=pipeArray.get(p);
				for(int j=0;j<no_commPipes;j++)
				{
					x[i*no_commPipes+j]=findPipe(i).getLength()*generator.nextDouble();
					System.out.println("Segment for pipe id: "+p.getPipeID()+" and segment no: "+j+" has initial length= "+x[i*no_commPipes+j]);
				}
			}
			// Starting point for flows.
			if (formulation == Formulation.PARALLEL_LINK) {
				for(Pipe p:this.inputPipes)
				{
					int i=pipeArray.get(p);
					double val=generator.nextDouble()*maxFlow*2-maxFlow;
					if(val>0)
					{
						x[i+no_pipes*no_commPipes] = val;
						x[i+no_pipes+no_pipes*no_commPipes]=0;
						orientation.put(findPipe(i),1);
					}
					else
					{
						x[i+no_pipes*no_commPipes] = 0;
						x[i+no_pipes+no_pipes*no_commPipes]=-val;
						orientation.put(findPipe(i),-1);
					}
					System.out.println("Flow for pipe id: "+p.getPipeID()+" has initial value= "+val);
				}
			}
			else if (formulation == Formulation.DISCRETE_SEGMENT) {
				for(Pipe p:this.inputPipes)
				{
					int i=pipeArray.get(p);
					x[i+no_pipes*no_commPipes]=generator.nextDouble()*maxFlow*2-maxFlow;
					orientation.put(findPipe(i),(int)Math.signum(x[i+no_pipes*no_commPipes]));
					System.out.println("Flow for pipe id: "+p.getPipeID()+" has initial value= "+x[i+no_pipes*no_commPipes]);
				}
			}
			else {
				for(Pipe p:orientation.keySet())
				{
					int i=pipeArray.get(p);
					if(orientation.get(p)==1)
						x[i+no_pipes*no_commPipes]=maxFlow*generator.nextDouble();
					else if(orientation.get(p)==-1)
						x[i+no_pipes*no_commPipes]=-maxFlow*generator.nextDouble();
					else
						x[i+no_pipes*no_commPipes]=maxFlow*2*generator.nextDouble()-maxFlow;

					System.out.println("Flow for pipe id: "+p.getPipeID()+" has initial value= "+x[i+no_pipes*no_commPipes]);
				}
			}
		}
		return true;
	}
	
	protected boolean eval_f(int n, double[] x, boolean new_x, double[] obj_value)
    {
//		System.out.println("Entering eval_f");
		assert n == this.n;
        obj_value[0]=0;
        for( int i = 0; i < no_pipes; i++ ) 
        {
        	for(int j=0;j<no_commPipes;j++)
        		obj_value[0]+=x[i*no_commPipes+j]*findCommPipe(j).getCost();
        }
        return true;
    }

	protected boolean eval_grad_f(int n, double[] x, boolean new_x, double[] grad_f)
    {
        assert n == this.n;
        
        // Gradient of objective wrt length segment variables.
        for(int i=0;i<no_pipes;i++)
        {
        	for(int j=0;j<no_commPipes;j++)
        		grad_f[i*no_commPipes+j]=findCommPipe(j).getCost();
        }
	
        // Gradient of objective wrt flow variables.
        if (formulation == Formulation.PARALLEL_LINK) {
	        for(int i=0;i<2*no_pipes;i++)
	     	   grad_f[i+no_commPipes*no_pipes]=0;
        }
        else {
        	for(int i=0;i<no_pipes;i++)
          	   grad_f[i+no_commPipes*no_pipes]=0;
        }
        return true;
    }
	
	
	protected boolean eval_g(int n, double[] x, boolean new_x, int m, double[] g)
    {
        assert n == this.n;
        assert m == this.m;
        int i,j,k;

        // Flow conservation constraints.
        for(i = 0; i <no_nodes; i++ ) 
	    {
        	g[i]=0;
        	for(j=0;j<no_pipes;j++) {
        		if (formulation == Formulation.PARALLEL_LINK) {
        			g[i] += flowDir[i][j]*(x[j+no_pipes*no_commPipes]-x[j+no_pipes+no_pipes*no_commPipes]);
        		}
        		else {
        			g[i] += flowDir[i][j]*(x[j+no_pipes*no_commPipes]);
        		}
        	}
		    g[i]-=findNode(i).getDemand();
		}
		  
        
		// Pressure drop across cycle constraints.
		for(i = 0; i <no_cycles; i++ ) 
		{
			g[i+no_nodes]=0;
			for(j=0;j<no_pipes;j++)
			{
				for(k=0;k<no_commPipes;k++)
				{
					if (formulation == Formulation.PARALLEL_LINK) {
						g[i+no_nodes] += cycleFlowDir[i][j]*10.68*x[j*no_commPipes+k]*(safe_pow(x[j+no_pipes*no_commPipes]*0.001,1.852)-safe_pow(x[j+no_pipes+no_pipes*no_commPipes]*0.001,1.852))/(safe_pow(findCommPipe(k).getRoughness(),1.852)*safe_pow(findCommPipe(k).getDiameter()*0.001, 4.87));
					}
					else {
						g[i+no_nodes] += cycleFlowDir[i][j]*10.68*x[j*no_commPipes+k]*Math.signum(x[j+no_pipes*no_commPipes])*Math.pow(Math.abs((x[j+no_pipes*no_commPipes]*0.001)/findCommPipe(k).getRoughness()),1.852)/Math.pow((findCommPipe(k).getDiameter()*0.001), 4.87);
					}
				}
			}
			
		}
		
	    // Headloss at a node constraints.
		for(i = 0; i < no_nodes; i++ ) 
		{
			g[i+no_nodes+no_cycles]=0;
			for(j=0;j<no_pipes;j++)
			{
				for(k=0;k<no_commPipes;k++)
				{
					if (formulation == Formulation.PARALLEL_LINK) {
						g[i+no_nodes+no_cycles] += primLinkFlow[i][j]*10.68*x[j*no_commPipes+k]*(safe_pow(x[j+no_pipes*no_commPipes]*0.001,1.852)-safe_pow(x[j+no_pipes+no_pipes*no_commPipes]*0.001,1.852))/(safe_pow(findCommPipe(k).getRoughness(),1.852)*safe_pow(findCommPipe(k).getDiameter()*0.001, 4.87));
					}
					else {
						g[i+no_nodes+no_cycles] += primLinkFlow[i][j]*10.68*x[j*no_commPipes+k]*Math.signum(x[j+no_pipes*no_commPipes])*Math.pow(Math.abs((x[j+no_pipes*no_commPipes]*0.001)/findCommPipe(k).getRoughness()),1.852)/Math.pow((findCommPipe(k).getDiameter()*0.001), 4.87);
					}
				}
			}
		}
		
		// Constraints on flow in link.
		for (i = 0; i <  no_pipes; i++) {
			if (formulation == Formulation.PARALLEL_LINK) {
				g[i+no_nodes*2+no_cycles] = x[i+no_pipes*no_commPipes]*x[i+no_pipes+no_pipes*no_commPipes];
			}
			else {
				g[i+no_nodes*2+no_cycles] = Math.abs(x[i+no_pipes*no_commPipes]);
			}
		}

		// Sum of link segments constraints.
		for(i=0;i<no_pipes;i++)
		{
			g[i+2*no_nodes+no_pipes+no_cycles]=0;
			for(j=0;j<no_commPipes;j++)
				g[i+2*no_nodes+no_pipes+no_cycles]+=x[i*no_commPipes+j];
		}
		return true;
    }
	
    protected boolean eval_jac_g(int n, double[] x, boolean new_x, int m, int nele_jac, int[] iRow, int[] jCol, double[] values)
    {
//    	System.out.println("Entering eval_jac_g");
    	if(values == null)
    	{
    		// Return the structure of the jacobian.
    		//Possible TODO: currently case dependant non-zero values are not being assigned
    		int i=0,j=0;//row and column counters
    		int nonZeroEle=0;//counter for the number of non zero entry 
    		
    		// Flow conservation constraints.
    		for(i=0;i<no_nodes;i++)
    		{
    			for(j=0;j<no_pipes;j++)
    			{
    				iRow[nonZeroEle]=i;
					jCol[nonZeroEle]=j+no_pipes*no_commPipes;
					nonZeroEle++;
    			}
    			if (formulation == Formulation.PARALLEL_LINK) {
	    			for(j=0;j<no_pipes;j++)
	    			{
	    				iRow[nonZeroEle]=i;
						jCol[nonZeroEle]=j+no_pipes+no_pipes*no_commPipes;
						nonZeroEle++;
	    			}
    			}
    		}
    		
    		// Pressure drop across cycle constraints.
    		for(i=0;i<no_cycles;i++)
    		{
    			if (formulation == Formulation.PARALLEL_LINK) {
					for(j=0;j<(no_commPipes+2)*no_pipes;j++)
					{
						iRow[nonZeroEle]=i+no_nodes;
						jCol[nonZeroEle]=j;
						nonZeroEle++;
					}
    			}
    			else {
    				for(j=0;j<(no_commPipes+1)*no_pipes;j++)
        			{
        				iRow[nonZeroEle]=i+no_nodes;
        				jCol[nonZeroEle]=j;
        				nonZeroEle++;
        			}
    			}
    		}
    		
    		// Headloss at node constraints.
    		for(i=0;i<no_nodes;i++)
    		{
    			if (formulation == Formulation.PARALLEL_LINK) {
	    			for(j=0;j<(no_commPipes+2)*no_pipes;j++)
	    			{
	    				iRow[nonZeroEle]=i+no_nodes+no_cycles;
	    				jCol[nonZeroEle]=j;
	    				nonZeroEle++;
	    			}
    			}
    			else {
    				for(j=0;j<(no_commPipes+1)*no_pipes;j++)
	    			{
	    				iRow[nonZeroEle]=i+no_nodes+no_cycles;
	    				jCol[nonZeroEle]=j;
	    				nonZeroEle++;
	    			}
    			}
    		}
    		
    		// Constraints on flow in link.
    		for(i=0;i<no_pipes;i++)
    		{
    			iRow[nonZeroEle]=i+no_nodes*2+no_cycles;
    			jCol[nonZeroEle]=i+no_pipes*no_commPipes;
    			nonZeroEle++;
    			if (formulation == Formulation.PARALLEL_LINK) {
	    			iRow[nonZeroEle]=i+no_nodes*2+no_cycles;
	    			jCol[nonZeroEle]=i+no_pipes+no_pipes*no_commPipes;
	    			nonZeroEle++;
    			}
    		}
    		
    		// Sum of link segments constraints.
    		for(i=0;i<no_pipes;i++)
    		{
    			for(j=0;j<no_commPipes;j++)
    			{
    				iRow[nonZeroEle]=i+no_nodes*2+no_cycles+no_pipes;
    				jCol[nonZeroEle]=i*no_commPipes+j;
    				nonZeroEle++;
    			}
    		}
    	}
    	else
    	{
    		// Return the values for the corresponding non zero elements in the jacobian.
    		int i=0,j=0;//row and column counters
    		int nonZeroEle=0;//counter for the number of non zero entry 
    		
    		// Flow conservation constraints.
    		for(i=0;i<no_nodes;i++)
    		{
    			for(j=0;j<no_pipes;j++)
    				values[nonZeroEle++]=flowDir[i][j];
    			if (formulation == Formulation.PARALLEL_LINK) {
	    			for(j=0;j<no_pipes;j++)
	    				values[nonZeroEle++]=-flowDir[i][j];
    			}
    		}
    		
    		// Pressure drop across cycle constraints.
    		for(i=0;i<no_cycles;i++)
    		{
    			for(j=0;j<no_pipes;j++)
    			{
    				for(int k=0;k<no_commPipes;k++)
    				{
    					if (formulation == Formulation.PARALLEL_LINK) {
    						values[nonZeroEle++]=cycleFlowDir[i][j]*10.68*(safe_pow(x[j+no_pipes*no_commPipes]*0.001,1.852)-safe_pow(x[j+no_pipes+no_pipes*no_commPipes]*0.001,1.852))/(safe_pow(findCommPipe(k).getRoughness(),1.852)*safe_pow(findCommPipe(k).getDiameter()*0.001, 4.87));
    					}
    					else {
    						values[nonZeroEle++]=cycleFlowDir[i][j]*10.68*Math.signum(x[j+no_pipes*no_commPipes])*Math.pow(Math.abs((x[j+no_pipes*no_commPipes]*0.001)/findCommPipe(k).getRoughness()),1.852)/Math.pow((findCommPipe(k).getDiameter()*0.001), 4.87);
    					}
    				}
    			}
    			if (formulation == Formulation.PARALLEL_LINK) {
	    			for(j=0;j<no_pipes;j++)
	    			{
	    				values[nonZeroEle]=0;
	    				for(int k=0;k<no_commPipes;k++)
						{
	    					values[nonZeroEle]+=cycleFlowDir[i][j]*10.68*x[j*no_commPipes+k]*1.852*0.001*safe_pow((x[j+no_pipes*no_commPipes]*0.001),0.852)/(safe_pow(findCommPipe(k).getRoughness(), 1.852)*safe_pow((findCommPipe(k).getDiameter()*0.001), 4.87));
						}
	    				nonZeroEle++;
	    			}
	    			for(j=0;j<no_pipes;j++)
	    			{
	    				values[nonZeroEle]=0;
	    				for(int k=0;k<no_commPipes;k++)
						{
	    					values[nonZeroEle]-=cycleFlowDir[i][j]*10.68*x[j*no_commPipes+k]*1.852*0.001*safe_pow((x[j+no_pipes+no_pipes*no_commPipes]*0.001),0.852)/(safe_pow(findCommPipe(k).getRoughness(), 1.852)*safe_pow((findCommPipe(k).getDiameter()*0.001), 4.87));
						}
	    				nonZeroEle++;
	    			}
    			}
    			else {
    				for(j=0;j<no_pipes;j++)
	    			{
	    				values[nonZeroEle]=0;
	    				for(int k=0;k<no_commPipes;k++)
	    				{
	    					values[nonZeroEle]+=cycleFlowDir[i][j]*10.68*x[j*no_commPipes+k]*1.852*0.001*Math.pow(Math.abs((x[j+no_pipes*no_commPipes]*0.001)),0.852)/(Math.pow(findCommPipe(k).getRoughness(), 1.852)*Math.pow((findCommPipe(k).getDiameter()*0.001), 4.87));
	    				}
	    				nonZeroEle++;
	    			}
    			}
    		}
    		
    		// Headloss at node constraints.
    		if (formulation == Formulation.PARALLEL_LINK) {
    			for(i=0;i<no_nodes;i++)
        		{
        			for(j=0;j<no_pipes;j++)
        			{
        				for(int k=0;k<no_commPipes;k++)
        				{
        					values[nonZeroEle++]=primLinkFlow[i][j]*10.68*(safe_pow(x[j+no_pipes*no_commPipes]*0.001,1.852)-safe_pow(x[j+no_pipes+no_pipes*no_commPipes]*0.001,1.852))/(safe_pow(findCommPipe(k).getRoughness(),1.852)*safe_pow(findCommPipe(k).getDiameter()*0.001, 4.87));
        				}
        			}
        			for(j=0;j<no_pipes;j++)
        			{
        				values[nonZeroEle]=0;
        				for(int k=0;k<no_commPipes;k++)
    					{
        					values[nonZeroEle]+=primLinkFlow[i][j]*10.68*x[j*no_commPipes+k]*1.852*0.001*safe_pow((x[j+no_pipes*no_commPipes]*0.001),0.852)/(safe_pow(findCommPipe(k).getRoughness(), 1.852)*safe_pow((findCommPipe(k).getDiameter()*0.001), 4.87));
    					}
        				nonZeroEle++;
        			}
        			for(j=0;j<no_pipes;j++)
        			{
        				values[nonZeroEle]=0;
        				for(int k=0;k<no_commPipes;k++)
    					{
        					values[nonZeroEle]-=primLinkFlow[i][j]*10.68*x[j*no_commPipes+k]*1.852*0.001*safe_pow((x[j+no_pipes+no_pipes*no_commPipes]*0.001),0.852)/(safe_pow(findCommPipe(k).getRoughness(), 1.852)*safe_pow((findCommPipe(k).getDiameter()*0.001), 4.87));
    					}
        				nonZeroEle++;
        			}
        		}
    		}
    		else {
    			for(i=0;i<no_nodes;i++)
        		{
        			for(j=0;j<no_pipes;j++)
        			{
        				for(int k=0;k<no_commPipes;k++)
        				{
        					values[nonZeroEle++]=primLinkFlow[i][j]*10.68*Math.signum(x[j+no_pipes*no_commPipes])*Math.pow(Math.abs((x[j+no_pipes*no_commPipes]*0.001)/findCommPipe(k).getRoughness()),1.852)/Math.pow((findCommPipe(k).getDiameter()*0.001), 4.87);
        				}
        			}
        			for(j=0;j<no_pipes;j++)
        			{
        				values[nonZeroEle]=0;
        				for(int k=0;k<no_commPipes;k++)
    					{
        					values[nonZeroEle]+=primLinkFlow[i][j]*10.68*x[j*no_commPipes+k]*1.852*0.001*Math.pow(Math.abs((x[j+no_pipes*no_commPipes]*0.001)),0.852)/(Math.pow(findCommPipe(k).getRoughness(), 1.852)*Math.pow((findCommPipe(k).getDiameter()*0.001), 4.87));
    					}
        				nonZeroEle++;
        			}
        		}
    		}
    		
    		
    		// Constraints on flow in link.
    		for(i=0;i<no_pipes;i++)
    		{
    			if (formulation == Formulation.PARALLEL_LINK) {
	    			values[nonZeroEle++]=x[i+no_pipes+no_pipes*no_commPipes];
	    			values[nonZeroEle++]=x[i+no_pipes*no_commPipes];
    			}
    			else {
    				values[nonZeroEle++]=Math.signum(x[i+no_pipes*no_commPipes]);
    			}
    		}
    		
    		// Sum of link segments constraints.
    		for(i=0;i<no_pipes;i++)
    		{
    			for(j=0;j<no_commPipes;j++)
    			{
    				values[nonZeroEle++]=1;
    			}
    		}
    	}
        
		return true;
    }
    protected boolean eval_h(int n, double[] x, boolean new_x, double obj_factor, int m, double[] lambda, boolean new_lambda, int nele_hess, int[] iRow, int[] jCol, double[] values)
    {
//    	System.out.println("Entering eval_h");
    	if(values==null)
    	{
    		// Return the structure of the hessian. Since it is symmetric return only the lower triangle.
    		int i,j;
    		int nonZeroEle=0;
    		
    		// Parallel link: differentiating with respect to length of the jth segment of the ith pipe followed by the flow in the ith pipe
    		// Discrete segment: differentiating with respect to the flow in the ith pipe followed by length of the jth segment of the ith pipe
    		for(i=0;i<no_pipes;i++)
    		{
    			for(j=0;j<no_commPipes;j++)
    			{
    				if (formulation == Formulation.PARALLEL_LINK) {
	    				iRow[nonZeroEle]=no_pipes*no_commPipes+i;
	    				jCol[nonZeroEle]=i*no_commPipes+j;
	    				nonZeroEle++;
	    				iRow[nonZeroEle]=no_pipes*no_commPipes+no_pipes+i;
	    				jCol[nonZeroEle]=i*no_commPipes+j;
    				}
    				else {
    					iRow[nonZeroEle]=i+no_pipes*no_commPipes;
    					jCol[nonZeroEle]=i*no_commPipes+j;
    				}
    				nonZeroEle++;
    			}
    		}
    		
    		// Parallel link: differentiating twice wrt flow 1 in ith pipe
    		// Discrete segment: differentiating twice with respect to flow in ith pipe
    		for(i=0;i<no_pipes;i++)
    		{
    			iRow[nonZeroEle]=i+no_pipes*no_commPipes;
    			jCol[nonZeroEle]=i+no_pipes*no_commPipes;
    			nonZeroEle++;
    		}
    		
    		if (formulation == Formulation.PARALLEL_LINK) {
	    		// Differentiating first wrt flow 1 in ith pipe, then with flow2 in this pipe.
	    		for(i=0;i<no_pipes;i++)
	    		{
	    			iRow[nonZeroEle]=i+no_pipes*no_commPipes+no_pipes;
	    			jCol[nonZeroEle]=i+no_pipes*no_commPipes;
	    			nonZeroEle++;
	    		}
	    		// Differentiating twice wrt flow 2 in ith pipe.
	    		for(i=0;i<no_pipes;i++)
	    		{
	    			iRow[nonZeroEle]=i+no_pipes*no_commPipes+no_pipes;
	    			jCol[nonZeroEle]=i+no_pipes*no_commPipes+no_pipes;
	    			nonZeroEle++;
	    		}
    		}
    	}
    	else
    	{
    		int i,j;
    		if (formulation == Formulation.PARALLEL_LINK) {
	    		for(i=0;i<no_pipes*(3+2*no_commPipes);i++)
	    		{
	    			values[i]=0;
	    		}
    		} else {
    			for(i=0;i<no_pipes*(1+no_commPipes);i++)
    			{
    				values[i]=0;
    			}
    		}
//    		System.out.println("Here");
    		
    		// Objective, and flow conservation constraints are linear, 0 hessian
    		
    		//////////////////////////  Pressure drop across cycle constraints.
    		// Parallel link: differentiating with respect to length of the jth segment of the ith pipe followed by the flow in the ith pipe
    		// Discrete segment: differentiating with respect to the flow in the ith pipe followed by length of the jth segment of the ith pipe
    		for(i=0;i<no_cycles;i++)
    		{
    			for(j=0;j<no_pipes;j++)
    			{
    				for(int k=0;k<no_commPipes;k++)
    				{
    					if (formulation == Formulation.PARALLEL_LINK) {
	    					//values[j*no_commPipes+k]+=lambda[i+no_nodes]*cycleFlowDir[i][j]*10.68*Math.signum(x[j+no_pipes*no_commPipes])*1.852*0.001*pow(Math.abs((x[j+no_pipes*no_commPipes]*0.001)),0.852)/(pow(findCommPipe(k).getRoughness(),1.852)*pow((findCommPipe(k).getDiameter()*0.001), 4.87));
	    					values[2*(j*no_commPipes+k)]+=lambda[i+no_nodes]*cycleFlowDir[i][j]*10.68*1.852*0.001*safe_pow(x[j+no_pipes*no_commPipes]*0.001,0.852)/(safe_pow(findCommPipe(k).getRoughness(),1.852)*safe_pow((findCommPipe(k).getDiameter()*0.001), 4.87));
	    					values[2*(j*no_commPipes+k)+1]-=lambda[i+no_nodes]*cycleFlowDir[i][j]*10.68*1.852*0.001*safe_pow(x[j+no_pipes+no_pipes*no_commPipes]*0.001,0.852)/(safe_pow(findCommPipe(k).getRoughness(),1.852)*safe_pow((findCommPipe(k).getDiameter()*0.001), 4.87));
    					}
    					else {
    						//values[j*no_commPipes+k]+=lambda[i+no_nodes]*cycleFlowDir[i][j]*10.68*Math.signum(x[j+no_pipes*no_commPipes])*1.852*0.001*Math.pow(Math.abs((x[j+no_pipes*no_commPipes]*0.001)),0.852)/(Math.pow(findCommPipe(k).getRoughness(),1.852)*Math.pow((findCommPipe(k).getDiameter()*0.001), 4.87));
    						values[j*no_commPipes+k]+=lambda[i+no_nodes]*cycleFlowDir[i][j]*10.68*1.852*0.001*Math.pow(Math.abs((x[j+no_pipes*no_commPipes]*0.001)),0.852)/(Math.pow(findCommPipe(k).getRoughness(),1.852)*Math.pow((findCommPipe(k).getDiameter()*0.001), 4.87));
    					}
    				}
    			}
    		}
//    		System.out.println("Here1");
    		
    		// Parallel link: differentiating twice with respect to flow 1 in ith pipe
    		// Discrete segment: differentiating twice with respect to flow in ith pipe
    		for(i=0;i<no_cycles;i++)
    		{
    			for(j=0;j<no_pipes;j++)
    			{
    				for(int k=0;k<no_commPipes;k++)
    				{
    					if (formulation == Formulation.PARALLEL_LINK) {
	    					values[j+2*no_pipes*no_commPipes]+=lambda[i+no_nodes]*cycleFlowDir[i][j]*10.68*x[j*no_commPipes+k]*1.852*0.852*0.001*0.001*safe_pow(x[j+no_pipes*no_commPipes]*0.001,-0.148)/(safe_pow(findCommPipe(k).getRoughness(),1.852)*safe_pow((findCommPipe(k).getDiameter()*0.001), 4.87));
	    					//values[j+no_pipes*no_commPipes]+=lambda[i+no_nodes]*cycleFlowDir[i][j]*10.68*x[j*no_commPipes+k]*1.852*0.852*0.001*0.001*pow(Math.abs((x[j+no_pipes*no_commPipes]*0.001)),-0.148)/(pow(findCommPipe(k).getRoughness(),1.852)*pow((findCommPipe(k).getDiameter()*0.001), 4.87));
    					}
    					else {
    						values[j+no_pipes*no_commPipes]+=lambda[i+no_nodes]*cycleFlowDir[i][j]*10.68*x[j*no_commPipes+k]*Math.signum(x[j+no_pipes*no_commPipes])*1.852*0.852*0.001*0.001*Math.pow(Math.abs((x[j+no_pipes*no_commPipes]*0.001)),-0.148)/(Math.pow(findCommPipe(k).getRoughness(),1.852)*Math.pow((findCommPipe(k).getDiameter()*0.001), 4.87));
    						//values[j+no_pipes*no_commPipes]+=lambda[i+no_nodes]*cycleFlowDir[i][j]*10.68*x[j*no_commPipes+k]*1.852*0.852*0.001*0.001*Math.pow(Math.abs((x[j+no_pipes*no_commPipes]*0.001)),-0.148)/(Math.pow(findCommPipe(k).getRoughness(),1.852)*Math.pow((findCommPipe(k).getDiameter()*0.001), 4.87));
    					}
    				}
    			}
    		}
//    		System.out.println("Here2");
    		
    		if (formulation == Formulation.PARALLEL_LINK) {
	    		// Differentiating twice wrt flow 2 in ith pipe.
	    		for(i=0;i<no_cycles;i++)
	    		{
	    			for(j=0;j<no_pipes;j++)
	    			{
	    				for(int k=0;k<no_commPipes;k++)
	    				{
	    					values[j+2*(no_pipes*no_commPipes+no_pipes)]-=lambda[i+no_nodes]*cycleFlowDir[i][j]*10.68*x[j*no_commPipes+k]*1.852*0.852*0.001*0.001*safe_pow(x[j+no_pipes*no_commPipes+no_pipes]*0.001,-0.148)/(safe_pow(findCommPipe(k).getRoughness(),1.852)*safe_pow((findCommPipe(k).getDiameter()*0.001), 4.87));
	    					//values[j+no_pipes*no_commPipes]+=lambda[i+no_nodes]*cycleFlowDir[i][j]*10.68*x[j*no_commPipes+k]*1.852*0.852*0.001*0.001*pow(Math.abs((x[j+no_pipes*no_commPipes]*0.001)),-0.148)/(pow(findCommPipe(k).getRoughness(),1.852)*pow((findCommPipe(k).getDiameter()*0.001), 4.87));
	    				}
	    			}
	    		}
    		}
//    		System.out.println("Here3");
    		
    		////////////////////////// Headloss at node constraints.
    		// Parallel link: differentiating with respect to length of the jth segment of the ith pipe followed by the flow in the ith pipe
    		// Discrete segment: differentiating with respect to the flow in the ith pipe followed by length of the jth segment of the ith pipe
    		for(i=0;i<no_nodes;i++)
    		{
    			for(j=0;j<no_pipes;j++)
    			{
    				for(int k=0;k<no_commPipes;k++)
    				{
    					if (formulation == Formulation.PARALLEL_LINK) {
	    					//values[j*no_commPipes+k]+=lambda[i+no_nodes]*cycleFlowDir[i][j]*10.68*Math.signum(x[j+no_pipes*no_commPipes])*1.852*0.001*pow(Math.abs((x[j+no_pipes*no_commPipes]*0.001)),0.852)/(pow(findCommPipe(k).getRoughness(),1.852)*pow((findCommPipe(k).getDiameter()*0.001), 4.87));
	    					values[2*(j*no_commPipes+k)]+=lambda[i+no_nodes+no_cycles]*primLinkFlow[i][j]*10.68*1.852*0.001*safe_pow(x[j+no_pipes*no_commPipes]*0.001,0.852)/(safe_pow(findCommPipe(k).getRoughness(),1.852)*safe_pow((findCommPipe(k).getDiameter()*0.001), 4.87));
	    					values[2*(j*no_commPipes+k)+1]-=lambda[i+no_nodes+no_cycles]*primLinkFlow[i][j]*10.68*1.852*0.001*safe_pow(x[j+no_pipes+no_pipes*no_commPipes]*0.001,0.852)/(safe_pow(findCommPipe(k).getRoughness(),1.852)*safe_pow((findCommPipe(k).getDiameter()*0.001), 4.87));
    					}
    					else {
    						//values[j*no_commPipes+k]+=lambda[i+no_nodes+no_cycles]*primLinkFlow[i][j]*10.68*Math.signum(x[j+no_pipes*no_commPipes])*1.852*0.001*Math.pow(Math.abs((x[j+no_pipes*no_commPipes]*0.001)),0.852)/(Math.pow(findCommPipe(k).getRoughness(),1.852)*Math.pow((findCommPipe(k).getDiameter()*0.001), 4.87));
    						values[j*no_commPipes+k]+=lambda[i+no_nodes+no_cycles]*primLinkFlow[i][j]*10.68*1.852*0.001*Math.pow(Math.abs((x[j+no_pipes*no_commPipes]*0.001)),0.852)/(Math.pow(findCommPipe(k).getRoughness(),1.852)*Math.pow((findCommPipe(k).getDiameter()*0.001), 4.87));
    					}
    				}
    			}
    		}
    		if (formulation == Formulation.PARALLEL_LINK) {
	    		// differentiating twice with respect to flow 1 in ith pipe
	    		for(i=0;i<no_nodes;i++)
	    		{
	    			for(j=0;j<no_pipes;j++)
	    			{
	    				for(int k=0;k<no_commPipes;k++)
	    				{
	    					values[j+2*no_pipes*no_commPipes]+=lambda[i+no_nodes+no_cycles]*primLinkFlow[i][j]*10.68*x[j*no_commPipes+k]*1.852*0.852*0.001*0.001*safe_pow(x[j+no_pipes*no_commPipes]*0.001,-0.148)/(safe_pow(findCommPipe(k).getRoughness(),1.852)*safe_pow((findCommPipe(k).getDiameter()*0.001), 4.87));
	    					//values[j+no_pipes*no_commPipes]+=lambda[i+no_nodes]*cycleFlowDir[i][j]*10.68*x[j*no_commPipes+k]*1.852*0.852*0.001*0.001*pow(Math.abs((x[j+no_pipes*no_commPipes]*0.001)),-0.148)/(pow(findCommPipe(k).getRoughness(),1.852)*pow((findCommPipe(k).getDiameter()*0.001), 4.87));
	    				}
	    			}
	    		}
    		}
//    		System.out.println("Here4");
    		
    		// Parallel link: differentiating twice wrt flow 2 in ith pipe
    		// Discrete segment: differentiating twice with respect to flow in ith pipe
    		for(i=0;i<no_nodes;i++)
    		{
    			for(j=0;j<no_pipes;j++)
    			{
    				for(int k=0;k<no_commPipes;k++)
    				{
    					if (formulation == Formulation.PARALLEL_LINK) {
	    					values[j+2*(no_pipes*no_commPipes+no_pipes)]-=lambda[i+no_nodes+no_cycles]*primLinkFlow[i][j]*10.68*x[j*no_commPipes+k]*1.852*0.852*0.001*0.001*safe_pow(x[j+no_pipes*no_commPipes+no_pipes]*0.001,-0.148)/(safe_pow(findCommPipe(k).getRoughness(),1.852)*safe_pow((findCommPipe(k).getDiameter()*0.001), 4.87));
	    					//values[j+no_pipes*no_commPipes]+=lambda[i+no_nodes]*cycleFlowDir[i][j]*10.68*x[j*no_commPipes+k]*1.852*0.852*0.001*0.001*pow(Math.abs((x[j+no_pipes*no_commPipes]*0.001)),-0.148)/(pow(findCommPipe(k).getRoughness(),1.852)*pow((findCommPipe(k).getDiameter()*0.001), 4.87));
    					}
    					else {
    						values[j+no_pipes*no_commPipes]+=lambda[i+no_nodes+no_cycles]*primLinkFlow[i][j]*10.68*x[j*no_commPipes+k]*Math.signum(x[j+no_pipes*no_commPipes])*1.852*0.852*0.001*0.001*Math.pow(Math.abs((x[j+no_pipes*no_commPipes]*0.001)),-0.148)/(Math.pow(findCommPipe(k).getRoughness(),1.852)*Math.pow((findCommPipe(k).getDiameter()*0.001), 4.87));
    					}
    				}
    			}
    		}
//    		System.out.println("Here5");
    		
    		if (formulation == Formulation.PARALLEL_LINK) {
    			////////////////////////// Constraints on flow in link.
	    		for(i=0;i<no_pipes;i++)
	    		{
	    			values[i+2*(no_pipes*no_commPipes)+no_pipes]+=lambda[i+2*no_nodes+no_cycles];
	    		}
	    		// Sum of link segments constraints are again linear, 0 hessian
    		}
    		else {
    			//Constraints on flow in link and sum of link segments constraints are again linear, 0 hessian
    		}
    		
    		
    	}
        return true;
    }
    
    private double safe_pow(double no,double exp)
    {
    	if(no<=0)
    		return 0;
    	else
    		return Math.pow(no, exp);
    }
    public int OptimizeNLP(double start)
    {
    	this.start=start;
    	return OptimizeNLP();
    }

    public int OptimizeNLP(long seed)
    {
    	this.seed=seed;
    	return this.OptimizeNLP();
    }
    public int OptimizeNLP(HashMap<Pipe,Integer> currOrient)
    {
//    	System.out.println("Called the orientation optimize function");
    	orientation=currOrient;
    	int ans=OptimizeNLP();
    	orientation=null;
    	return ans;
    }
}