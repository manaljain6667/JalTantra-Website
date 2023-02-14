package optimizer;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Set;

public class redGraph {
	public HashMap<Integer,ArrayList<edgeGroup>> adjList;
	/**source for the original graph*/
	int source;
	final Boolean parentIsNetwork;
	redGraph parentGraph;
	redGraph simpleGraph;
	Set<edgeGroup> bridges;
	ArrayList<Integer> reducedSources;
	/**list of all sources of each connected component of the graph*/
	ArrayList<Integer> componentSources;///make a set?
	HashMap<edgeGroup,Integer> bridgeOrient;
	Boolean end=true;
	int maxOrientAtEachStep=10;
	/**redGraph equivalent of given graph*/
	redGraph(int source,ArrayList<Integer> componentSources,int maxOrient)
	{
		this.source=source;
		this.parentIsNetwork=true;
		this.componentSources=componentSources;
		parentGraph=null;
		bridges=null;
		simpleGraph=null;
		reducedSources=null;
		adjList=new HashMap<Integer,ArrayList<edgeGroup>>();
		bridgeOrient=new HashMap<edgeGroup,Integer>();
		maxOrientAtEachStep=maxOrient;
	}
	/**if reduced version of existing redGraph, need to have reference to the original graph*/
	redGraph(int source,redGraph parentGraph,ArrayList<Integer> componentSources,int maxOrient)
	{
		this.source=source;
		this.parentIsNetwork=false;
		this.parentGraph=parentGraph;
		this.componentSources=componentSources;
		bridges=null;
		simpleGraph=null;
		reducedSources=null;
		adjList=new HashMap<Integer,ArrayList<edgeGroup>>();
		bridgeOrient=new HashMap<edgeGroup,Integer>();
		maxOrientAtEachStep=maxOrient;
	}
	
	
	/**dfs for finding bridges by computing low, level values for each node*/
	private void graphDFS(HashMap<Integer,Integer> low,HashMap<Integer,Integer> level, int this_level,int curr_node,HashMap<Integer,Boolean> visit,Set<edgeGroup> bridges)
	{
		int nextNode;
		visit.put(curr_node,true);
		level.put(curr_node,this_level);
		low.put(curr_node,this_level);
		for(edgeGroup p:simpleGraph.adjList.get(curr_node))
		{
//			if(p.getStartNode()==curr_node)
//				nextNode=p.getEndNode();
//			else
//				nextNode=p.getStartNode();
			nextNode=p.otherNode(curr_node);
			
			if(!visit.get(nextNode))
			{
				graphDFS(low,level,this_level+1,nextNode,visit,bridges);
				low.put(curr_node,Math.min(low.get(nextNode), low.get(curr_node)));
				if(low.get(nextNode)>this_level)
				{//the edge is a bridge
					bridges.add(p);
				}
			}
			else if(1+level.get(nextNode)!=level.get(curr_node))
			{//back edge
				low.put(curr_node, Math.min(level.get(nextNode), low.get(curr_node)));
			}
			else
			{//tree edge from child to parent
				
			}
			
			
		}
	}
	
	/**all bridges in all connected components. 
	 * do nothing if they are already made*/
	private void bridgeSet()
	{
		if(bridges!=null)
			return;
		HashMap<Integer,Integer> low=new HashMap<Integer,Integer>(),level=new HashMap<Integer,Integer>();
		HashMap<Integer,Boolean> visit=new HashMap<Integer,Boolean>();
		bridges=new HashSet<edgeGroup>();
		reducedSources=new ArrayList<Integer>(componentSources);
		for(Integer i: simpleGraph.adjList.keySet())
		{
			visit.put(i, false);
		}
		for(int n:componentSources)
		{
			if(!visit.get(n))
			{
				graphDFS(low,level,0,n,visit,bridges);//level initialised to 0 for each connected component
			}
		}
		
		//bridges made
		
		for(edgeGroup b:bridges)
		{
			if(level.get(b.getNode1())>level.get(b.getNode2()))
			{//node1 is further from source so it will be in the new conn comp and the source for that
				bridgeOrient.put(b, b.towards(b.getNode1()));
				reducedSources.add(b.getNode1());
			}
			else
			{
				bridgeOrient.put(b, b.towards(b.getNode2()));
				reducedSources.add(b.getNode2());
			}
		}
		return;
	}
	
	/**give a redGraph. 
	 * first make the simpleGraph, then check if it is reducible.
	 * If it is, then reduce and return the reduction of the reduced graph recursively*/
	public redGraph reduce()
	{
		if(!reducible())
			return this;
		toSimple();
		bridgeSet();
		redGraph postReduction=new redGraph(source,this,reducedSources,maxOrientAtEachStep);
		HashMap<Integer,Boolean> visit=new HashMap<Integer,Boolean>();
		HashMap<Integer,Integer> level=new HashMap<Integer,Integer>();
		Set<edgeGroup> remEdge=new HashSet<edgeGroup>();
		for(int n:simpleGraph.adjList.keySet())
		{
			ArrayList<edgeGroup> newarli=new ArrayList<edgeGroup>();
			postReduction.adjList.put(n,newarli);
			visit.put(n,false);
			level.put(n,-1);
		}
		int highestEdge=-1;
		for(int s:componentSources)
		{
			level.put(s, 0);
//			edgeGroup feg=new edgeGroup(highestEdge+1,source);
//			highestEdge+=1;
			highestEdge=redGraphDFS(postReduction,null,s,visit,false,highestEdge,level,remEdge,bridges);
		}
		
		
		
		return postReduction.reduce();
	}
	
	/**makes a simpleGraph out of the current one by collapsing multiple edges and removing loops. 
	 * if it is already made, then do nothing*/
	private void toSimple()
	{
		//to check if the node is already added
//		HashMap<Integer,Set<Integer>> simpleAdjList=new HashMap<Integer,Set<Integer>>();
		if(simpleGraph!=null)
			return;
		int highestID=0;
		simpleGraph=new redGraph(source,componentSources,maxOrientAtEachStep);
		HashMap<Integer,HashMap<Integer,edgeGroup>> corresEdge=new HashMap<Integer,HashMap<Integer,edgeGroup>>();
		for(int i:adjList.keySet())
		{
//			Set<Integer> newset=new HashSet<Integer>();
			ArrayList<edgeGroup> newarli=new ArrayList<edgeGroup>();
//			simpleAdjList.put(i,newset);
			simpleGraph.adjList.put(i, newarli);
			HashMap<Integer,edgeGroup> bla=new HashMap<Integer,edgeGroup>();
			corresEdge.put(i, bla);
		}
		for(int i:adjList.keySet())
		{
			//the edge in simple graph which connects i to any neighbour
			
			for(edgeGroup p:adjList.get(i))
			{
				int otherEnd=p.otherNode(i);
				if(otherEnd!=i)
				{//no self loop
					if(!corresEdge.get(i).containsKey(otherEnd))
					{//introduce a new edge
//						simpleAdjList.get(i).add(otherEnd);
//						simpleAdjList.get(otherEnd).add(i);
						edgeGroup neweg=new edgeGroup(highestID+1,i,otherEnd);
						neweg.addEdge(p);
						highestID+=1;
						simpleGraph.adjList.get(i).add(neweg);
						simpleGraph.adjList.get(otherEnd).add(neweg);
						corresEdge.get(i).put(otherEnd, neweg);
						corresEdge.get(otherEnd).put(i, neweg);
					}
					else
					{//add p to existing edge
						corresEdge.get(i).get(otherEnd).addEdge(p);
					}
					
					
				}
			}
		}
		
		return;
	}
//	private void makeSimpleGraph()
//	{
//		simpleGraph=new redGraph(source);
//		for(int i:adjList.keySet())
//		{
//			ArrayList<edgeGroup> newarli=new ArrayList<edgeGroup>();
//			simpleGraph.adjList.put(i, newarli);
//			for(edgeGroup p:)
//		}
//	}
	
	/**return simple graph with added nodes to represent ignored nodes
	 * if n1-(l)-n2 is such that l represents an edges that replaces some nodes,
	 * this returns n1-n'-n2 that is, adds a node in between to allow n' to be a sink*/
	private HashMap<Integer,Set<Integer>> toDerivedSimple(HashMap<Integer,Integer> extraNode)
	{
		HashMap<Integer,Set<Integer>> simpleAdjList=new HashMap<Integer,Set<Integer>>();
		//rem node is list of nodes whose connection will be broken and new nodes introduced
		ArrayList<Integer> remNode1=new ArrayList<Integer>(),remNode2=new ArrayList<Integer>(),remPipe=new ArrayList<Integer>();
		int highestNode=-1;
		for(int n:adjList.keySet())
			if(n>highestNode)
				highestNode=n;//highest node in the simpleGraph
//		highestNode+=1;
		for(int i:adjList.keySet())
		{
			Set<Integer> newarli=new HashSet<Integer>();
			simpleAdjList.put(i,newarli);
			for(edgeGroup p:adjList.get(i))
			{
				int otherEnd=p.otherNode(i);
				if(otherEnd!=i)
				{//if this is a self loop then ignore it, take care of self loops later as their orientations are fixed to be 0
					simpleAdjList.get(i).add(otherEnd);
					if(p.getMappingSize()>1)
					{//atleast 1 node was ignored
						remNode1.add(i);
						remNode2.add(otherEnd);
						remPipe.add(p.getID());
					}
				}
				
			}
		}
		int l=remNode1.size();
		for(int i=0;i<l;i++)
		{
			if(simpleAdjList.get(remNode1.get(i)).contains(remNode2.get(i)))
			{
				simpleAdjList.get(remNode1.get(i)).remove(remNode2.get(i));
				simpleAdjList.get(remNode2.get(i)).remove(remNode1.get(i));
				highestNode+=1;
				extraNode.put(remPipe.get(i),highestNode);
				simpleAdjList.get(remNode1.get(i)).add(highestNode);
				simpleAdjList.get(remNode2.get(i)).add(highestNode);
				Set<Integer> newarli=new HashSet<Integer>();
				newarli.add(remNode1.get(i));
				newarli.add(remNode2.get(i));
				simpleAdjList.put(highestNode,newarli);
			}
		}
		return simpleAdjList;
	}
	
	/**return all orientations of the original graph by going recursively. 
	 * here we are given all orientations for the graph obtained by reducing this*/
	private ArrayList<HashMap<edgeGroup,Integer>> enumerateOrientations(ArrayList<HashMap<edgeGroup,Integer>> childOrientations)
	{
		ArrayList<HashMap<edgeGroup,Integer>> simpOrient=fromChildOrient(childOrientations);
		ArrayList<HashMap<edgeGroup,Integer>> multOrient=simpToMultOrientations(simpOrient);
		if(parentIsNetwork)
			return multOrient;
		else
		{
			return parentGraph.enumerateOrientations(multOrient);
		}
	}
	
	/**return all orientations of the original graph by going recursively*/
	public ArrayList<HashMap<edgeGroup,Integer>> enumerateOrientations()
	{
		ArrayList<HashMap<edgeGroup,Integer>> simpOrient=makeSimpleOrientations();
		ArrayList<HashMap<edgeGroup,Integer>> multOrient=simpToMultOrientations(simpOrient);
		if(parentIsNetwork)//go up no more
			return multOrient;
		else
		{
			return parentGraph.enumerateOrientations(multOrient);
		}
	}
	
	/**given all orientations of the reduced multiple graph, make all orientations of the simple graph*/
	private ArrayList<HashMap<edgeGroup,Integer>> fromChildOrient(ArrayList<HashMap<edgeGroup,Integer>> childOrientations)
	{
		ArrayList<HashMap<edgeGroup,Integer>> orientations=new ArrayList<HashMap<edgeGroup,Integer>>();
		for(HashMap<edgeGroup,Integer> childOrient:childOrientations)
		{
			//all orientations possible from the given orientation
			ArrayList<HashMap<edgeGroup,Integer>> listNewOrient=new ArrayList<HashMap<edgeGroup,Integer>>();
			//for each edge group, list of possible orientations
			ArrayList<ArrayList<HashMap<edgeGroup,Integer>>> lolOrient=new ArrayList<ArrayList<HashMap<edgeGroup,Integer>>>();
			for(edgeGroup p: childOrient.keySet())
			{
				ArrayList<HashMap<edgeGroup,Integer>> egOrient=new ArrayList<HashMap<edgeGroup,Integer>>();
				//edges(among those mapped by p) connected to each node
				HashMap<Integer,Set<edgeGroup>> intNodes=new HashMap<Integer,Set<edgeGroup>>();
				for(edgeGroup p1:p.getMapping())
				{
					int node1=p1.getNode1(),node2=p1.getNode2();
					if(intNodes.containsKey(node1))
					{
						intNodes.get(node1).add(p1);
					}
					else
					{
						Set<edgeGroup> negs=new HashSet<edgeGroup>();
						negs.add(p1);
						intNodes.put(node1, negs);
						
					}
					if(intNodes.containsKey(node2))
					{
						intNodes.get(node2).add(p1);
					}
					else
					{
						Set<edgeGroup> negs=new HashSet<edgeGroup>();
						negs.add(p1);
						intNodes.put(node2, negs);
						
					}
				}
				if(childOrient.get(p)==1)
				{
					if(p.getNode1()!=p.getNode2())
						egOrient.add(orientSink(p.getNode2(),p.getMapping(),intNodes,-1));
					else
						egOrient.add(orientSink(p.getNode2(),p.getMapping(),intNodes,p.getNode1()));
				}
					
				else if(childOrient.get(p)==-1)
				{
					if(p.getNode1()!=p.getNode2())
						egOrient.add(orientSink(p.getNode1(),p.getMapping(),intNodes,-1));
					else
						egOrient.add(orientSink(p.getNode1(),p.getMapping(),intNodes,p.getNode1()));
				}
				else
				{//every node except the end ones are potential sinks
					for(int n:intNodes.keySet())
					{
						if(p.getNode1()!=p.getNode2())
						{
							if(intNodes.get(n).size()>1)
							{//not an end node
								egOrient.add(orientSink(n,p.getMapping(),intNodes,-1));
							}
						}
						else
						{
							if(n!=p.getNode1())
							{//not an end node
								egOrient.add(orientSink(n,p.getMapping(),intNodes,p.getNode1()));
							}
						}
					}
					//every component edge made of two or more edges can also have a sink
					for(edgeGroup simpEg:p.getMapping())
					{
						Boolean canBe0=false;
						for(edgeGroup compEg:simpEg.getMapping())
						{
							if(compEg.getMappingSize()>=2)
							{
								canBe0=true;
								break;
							}
						}
						if(canBe0)
						{
							if(p.getNode1()==p.getNode2())
								egOrient.add(orientEdgeSink(simpEg,p.getMapping(),intNodes,p.getNode1()));
							else
								egOrient.add(orientEdgeSink(simpEg,p.getMapping(),intNodes,-1));
						}
							
					}
				}
				lolOrient.add(egOrient);
				
			}
			for(edgeGroup br:bridgeOrient.keySet())
			{
				ArrayList<HashMap<edgeGroup,Integer>> bror=new ArrayList<HashMap<edgeGroup,Integer>>();
				HashMap<edgeGroup,Integer> bro=new HashMap<edgeGroup,Integer>();
				bro.put(br, bridgeOrient.get(br));
				bror.add(bro);
				lolOrient.add(bror);
			}
			listNewOrient=cartProd(lolOrient,maxOrientAtEachStep);
			orientations.addAll(listNewOrient);
		}
		System.out.println("In fromChildOrient(), generated "+orientations.size()+" number of orientations");
		orientations=reduceSize(orientations);
		System.out.println("Return "+orientations.size()+" number of orientations");
		return orientations;
	}
	
	/**return orientation setting sink to be in sinkEdge among edgeSet*/
	private HashMap<edgeGroup,Integer> orientEdgeSink(edgeGroup sinkEdge,Set<edgeGroup> edgeSet,HashMap<Integer,Set<edgeGroup>> intNodes,int loopVert)
	{
		HashMap<edgeGroup,Integer> orientation=new HashMap<edgeGroup,Integer>();
		edgeGroup currEdge;
		int curr=sinkEdge.getNode1(),other;
		Boolean done=false;
		Iterator<edgeGroup> it=intNodes.get(curr).iterator();
		currEdge=it.next();
		if(currEdge==sinkEdge)
		{
			if(curr!=loopVert && intNodes.get(curr).size()>1)
			{//go forward only in that case
				currEdge=it.next();
			}
			else
				done=true;
		}
		if(!done)
		{
			other=currEdge.otherNode(curr);
			orientation.put(currEdge, currEdge.towards(curr));
			while(other!=loopVert && intNodes.get(other).size()>1)
			{
				curr=other;
				it=intNodes.get(other).iterator();
				edgeGroup neweg=it.next();
				if(neweg==currEdge)
					neweg=it.next();
				currEdge=neweg;
				other=currEdge.otherNode(curr);
				orientation.put(currEdge, currEdge.towards(curr));
			}
		}
		
		
		curr=sinkEdge.getNode2();
		it=intNodes.get(curr).iterator();
		currEdge=it.next();
		done=false;
		if(currEdge==sinkEdge)
		{
			if(curr!=loopVert && intNodes.get(curr).size()>1)
			{//go forward only in that case
				currEdge=it.next();
			}
			else
				done=true;
		}
		if(!done)
		{
			other=currEdge.otherNode(curr);
			orientation.put(currEdge, currEdge.towards(curr));
			while(other!=loopVert && intNodes.get(other).size()>1)
			{
				curr=other;
				it=intNodes.get(other).iterator();
				edgeGroup neweg=it.next();
				if(neweg==currEdge)
					neweg=it.next();
				currEdge=neweg;
				other=currEdge.otherNode(curr);
				orientation.put(currEdge, currEdge.towards(curr));
			}
		}
		
		orientation.put(sinkEdge, 0);
		
		return orientation;
	}
	
	
	/**set sink as the sink of edges in edgeSet. Make an orientation accordingly*/
	private HashMap<edgeGroup,Integer> orientSink(int sink,Set<edgeGroup> edgeSet,HashMap<Integer,Set<edgeGroup>> intNodes,int loopVert)
	{
		HashMap<edgeGroup,Integer> orientation=new HashMap<edgeGroup,Integer>();
//		HashMap<Integer,Set<edgeGroup>> _intNodes=new HashMap<Integer,Set<edgeGroup>>(intNodes);
		for(edgeGroup p:intNodes.get(sink))
		{
			int other=p.otherNode(sink),curr=sink;
			edgeGroup curreg=p;
			orientation.put(p, p.towards(sink));
			if(loopVert!=-1)
			{
				while(other!=loopVert)
				{
					curr=other;
					Iterator<edgeGroup> it=intNodes.get(other).iterator();
					edgeGroup neweg=it.next();
					if(neweg==curreg)
						neweg=it.next();
					curreg=neweg;
					other=curreg.otherNode(curr);
					//all edge directed towards the sink
					orientation.put(curreg,curreg.towards(curr));
				}
			}
			else
			{
				while(intNodes.get(other).size()>1)
				{
					curr=other;
					Iterator<edgeGroup> it=intNodes.get(other).iterator();
					edgeGroup neweg=it.next();
					if(neweg==curreg)
						neweg=it.next();
					curreg=neweg;
					other=curreg.otherNode(curr);
					//all edge directed towards the sink
					orientation.put(curreg,curreg.towards(curr));
				}
			}
			
			
		}
		return orientation;
	}
	
	
	/**given orientations of simple graph, make multiple graph's orientations*/
	private ArrayList<HashMap<edgeGroup,Integer>> simpToMultOrientations(ArrayList<HashMap<edgeGroup,Integer>> simpOrientations)
	{
//		ArrayList<HashMap<edgeGroup,Integer>> simpOrientations=makeSimpleOrientations();
		ArrayList<HashMap<edgeGroup,Integer>> orientations=new ArrayList<HashMap<edgeGroup,Integer>>();
		ArrayList<edgeGroup> loops=new ArrayList<edgeGroup>();
		for(int n:adjList.keySet())
		{
			for(edgeGroup p:adjList.get(n))
			{
				if(p.otherNode(n)==n)
					loops.add(p);
			}
		}
		for(HashMap<edgeGroup,Integer> simpOrient:simpOrientations)
		{
			//list of list of orientations for each simple graph's edge
			ArrayList<ArrayList<HashMap<edgeGroup,Integer>>> resOrientations=new ArrayList<ArrayList<HashMap<edgeGroup,Integer>>>();
			Boolean throwOrient=false;
			for(edgeGroup p: simpOrient.keySet())
			{
				ArrayList<HashMap<edgeGroup,Integer>> egOrient=new ArrayList<HashMap<edgeGroup,Integer>>();
				egOrient=edgeMultOrientations(p.getMapping().iterator(),p.orientToNode(simpOrient.get(p)),egOrient);
				if(egOrient==null)
				{
					throwOrient=true;
					break;
				}
				resOrientations.add(egOrient);
			}
			if(throwOrient)
				continue;
			for(edgeGroup p:loops)
			{//make a hashmap and add it to a list and add that to resOrientations
				HashMap<edgeGroup,Integer> loop=new HashMap<edgeGroup,Integer>();
				loop.put(p, 0);
				ArrayList<HashMap<edgeGroup,Integer>> loopOr=new ArrayList<HashMap<edgeGroup,Integer>>();
				loopOr.add(loop);
				resOrientations.add(loopOr);
			}
			//cartesian product of the lists
			ArrayList<HashMap<edgeGroup,Integer>> blankCP=new ArrayList<HashMap<edgeGroup,Integer>>();
			
			//this is all orientations which can be generated from simpOrient orientation of simple graph
			ArrayList<HashMap<edgeGroup,Integer>> modifiedOrientations=cartProd(resOrientations,maxOrientAtEachStep);
			orientations.addAll(modifiedOrientations);
		}
//		for(HashMap<edgeGroup,Integer> someOrient:orientations)//since all loops were removed in the simple graph, corresponding orientations are introduced here
//			someOrient.putAll(loopOrient);
		System.out.println("In simpToMultOrientations(), generated "+orientations.size()+" number of orientations");
		orientations=reduceSize(orientations);
		System.out.println("Return "+orientations.size()+" number of orientations");
		return orientations;
	}
	
	private ArrayList<HashMap<edgeGroup,Integer>> reduceSize(ArrayList<HashMap<edgeGroup,Integer>> ip)
	{
		ArrayList<HashMap<edgeGroup,Integer>> op=new ArrayList<HashMap<edgeGroup,Integer>>();
		
		//only look at a total of maxOrientAtEachStep generated orientations
		int skip=ip.size()/maxOrientAtEachStep;
		skip=Math.max(skip, 1);
		for(int i=0;i<ip.size();i+=skip)
			op.add(ip.get(i));
		return op;
	}
	
	private long getTotalOrient(ArrayList<ArrayList<HashMap<edgeGroup,Integer>>> lists)
	{
		long prod=1;
		for(int i=0;i<lists.size();i++)
			prod*=lists.get(i).size();
		return prod;
//		int n=lists.size();
//		ArrayList<Double> revProd=new ArrayList<Double>();
//		int ind=n-1;
//		revProd.add((double)lists.get(ind).size());
//		ind--;
//		for(;ind>=0;ind--)
//		{
//			revProd.add(revProd.get(n-2-ind)*lists.get(ind).size());
//		}
//		ArrayList<Double> sizeProd=new ArrayList<Double>();
//		for(int i=n-1;i>=0;i--)
//		{
//			sizeProd.add(revProd.get(i));
//		}
//		return sizeProd;
	}
	
	private ArrayList<Integer> getSizeList(ArrayList<ArrayList<HashMap<edgeGroup,Integer>>> lists)
	{
		ArrayList<Integer> sizeList=new ArrayList<Integer>();
		for(int i=0;i<lists.size();i++)
			sizeList.add(lists.get(i).size());
		return sizeList;
	}
	
	/* for given orientationNo, what orientation in each list will be selected?*/
	private ArrayList<Integer> getIndexList(long orientationNo,ArrayList<Integer> sizeList)
	{
		ArrayList<Integer> indList=new ArrayList<Integer>();
		for(int i=0;i<sizeList.size();i++)
		{
			indList.add((int)orientationNo%sizeList.get(i));
			orientationNo=orientationNo/sizeList.get(i);
		}
		return indList;
	}
	
	/**given list of list of map, choose one map from each list and merge them*/
	private ArrayList<HashMap<edgeGroup,Integer>> cartProd(ArrayList<ArrayList<HashMap<edgeGroup,Integer>>> lists,int maxAllowed)
	{
		System.out.println("Latest version of cartProd being used");
		// total number of possible orientations that will result
		long totSize=getTotalOrient(lists);
		ArrayList<Integer> sizeList=getSizeList(lists);
		
		//only looks for around maxAllowed number of possibilities
		long skip=totSize/maxAllowed;
		skip=Math.max(skip, 1);
		ArrayList<HashMap<edgeGroup,Integer>> combinedOrientations=new ArrayList<HashMap<edgeGroup,Integer>>();
		long orientationNo=0;
		for(;orientationNo<totSize;orientationNo+=skip)
		{
			ArrayList<Integer> indexList=getIndexList(orientationNo,sizeList);
			HashMap<edgeGroup,Integer> combined=new HashMap<edgeGroup,Integer>();
			for(int i=0;i<indexList.size();i++)
				combined.putAll(lists.get(i).get(indexList.get(i)));
			combinedOrientations.add(combined);
		}
		return combinedOrientations;
		
//		if(listNo==lists.size())
//			return cpSoFar;
//		
//		ArrayList<HashMap<edgeGroup,Integer>> newCP=new ArrayList<HashMap<edgeGroup,Integer>>();
//		if(cpSoFar.size()==0)
//		{
//			for(HashMap<edgeGroup,Integer> entry:lists.get(listNo))
//			{
//				newCP.add(entry);
//			}
//		}
//		for(HashMap<edgeGroup,Integer> cp:cpSoFar)
//		{
//			for(HashMap<edgeGroup,Integer> entry:lists.get(listNo))
//			{
//				HashMap<edgeGroup,Integer> newEntry=new HashMap<edgeGroup,Integer>(cp);
//				newEntry.putAll(entry);
//				newCP.add(newEntry);
//			}
//		}
//		
//		
//		return cartProd(listNo+1,lists,newCP);
	}
	
	/**if orientation of an edge of simple graph is +1 or -1, 
	 * then the multiple graph's corresponding edges can have orientations +1 or 0, except for all having 0.
	 *  this returns list of all possibilities of multiple graph's orientations. 
	 * We are given iterator for the set of edges of the multiple graph,orientation for the simple graph's edge*/
	private ArrayList<HashMap<edgeGroup,Integer>> edgeMultOrientations(Iterator<edgeGroup> memberIt,int towNode,ArrayList<HashMap<edgeGroup,Integer>> egOrient)
	{
		edgeGroup thiseg;
		if(memberIt.hasNext())//memberIt has moved to the next member
			thiseg=(edgeGroup)memberIt.next();
		else
			return egOrient;
		
		ArrayList<HashMap<edgeGroup,Integer>> newEgOrient=new ArrayList<HashMap<edgeGroup,Integer>>();
		if(egOrient.size()==0)
		{
			HashMap<edgeGroup,Integer> newOrient=new HashMap<edgeGroup,Integer>();
			if(towNode==-1)
			{
				if(thiseg.getMappingSize()<=1)
					return null;
				newOrient.put(thiseg, 0);
				newEgOrient.add(newOrient);
			}
			else
			{
				
				if(memberIt.hasNext() && thiseg.getMappingSize()>1)
				{
					newOrient.put(thiseg, 0);
					newEgOrient.add(newOrient);
				}
				newOrient=new HashMap<edgeGroup,Integer>();
				newOrient.put(thiseg, thiseg.towards(towNode));
				newEgOrient.add(newOrient);
			}
		}
		for(HashMap<edgeGroup,Integer> oldOrient:egOrient)
		{
			HashMap<edgeGroup,Integer> newOrient=new HashMap<edgeGroup,Integer>(oldOrient);
			if(towNode==-1)
			{
				if(thiseg.getMappingSize()<=1)
					return null;
				newOrient.put(thiseg, 0);
				newEgOrient.add(newOrient);
			}
			else
			{
				newOrient.put(thiseg, 0);
				if(thiseg.getMappingSize()>1)
				{
					if(!memberIt.hasNext())
					{
						Boolean all0=true;
						for(edgeGroup po:oldOrient.keySet())
							if(oldOrient.get(po)!=0)
							{
								all0=false;
								break;
							}
						if(!all0)//all cant be 0
							newEgOrient.add(newOrient);
					}
					else
						newEgOrient.add(newOrient);
				}
				
				newOrient=new HashMap<edgeGroup,Integer>(oldOrient);
				newOrient.put(thiseg, thiseg.towards(towNode));
				newEgOrient.add(newOrient);
			}
		}
		return edgeMultOrientations(memberIt,towNode,newEgOrient);
	}
	
	/**orientations of the simple graph. 
	 * first make a graph that accounts for ignored nodes, 
	 * then group nodes according to connected component and 
	 * make all permutations of individual groups to determine flows in that graph. 
	 * Finally use those flows to determine flows in the simple graph corresponding to current graph*/
	private ArrayList<HashMap<edgeGroup,Integer>> makeSimpleOrientations()
	{
		//list of orientations. each orientation is a list of nodes, in topological sorted order
		ArrayList<ArrayList<Integer>> simpOrientations=new ArrayList<ArrayList<Integer>>();
		ArrayList<HashMap<edgeGroup,Integer>> orientations=new ArrayList<HashMap<edgeGroup,Integer>>();
		HashMap<Integer,Integer> extraNode=new HashMap<Integer,Integer>();
		HashMap<Integer,Set<Integer>> derSimpGraph=simpPermute(simpOrientations,extraNode);
		for(ArrayList<Integer> simpOrient:simpOrientations)
		{
			HashMap<Integer,Integer> revSimp=revOrient(simpOrient);
			HashMap<edgeGroup,Integer> orient=new HashMap<edgeGroup,Integer>();
			//for every node, orientation for every other neighbour node
			HashMap<Integer,HashMap<Integer,Integer>> orientMap=new HashMap<Integer,HashMap<Integer,Integer>>();
			for(int n:derSimpGraph.keySet())
			{
				if(simpleGraph.adjList.containsKey(n))
				{
					orientMap.put(n,new HashMap<Integer,Integer>());
					for(int o:derSimpGraph.get(n))
					{
						if(simpleGraph.adjList.containsKey(o))
						{
							if(revSimp.get(o)>revSimp.get(n))
								orientMap.get(n).put(o,+1);//out from n
							else
								orientMap.get(n).put(o,-1);
						}
						else
						{// this edge was created to allow for ignored nodes
							int ot=o;
							for(int j:derSimpGraph.get(o))
								if(j!=n)
								{
									ot=j;
									break;
								}
							if(revSimp.get(o)>revSimp.get(ot) && revSimp.get(o)>revSimp.get(n))
								orientMap.get(n).put(ot,0);
							else if(revSimp.get(n)>revSimp.get(ot))
								orientMap.get(n).put(ot,-1);
							else //(revSimp.get(n)<revSimp.get(ot))
								orientMap.get(n).put(ot,1);
							
						}
					}
				}
			}
			for(int n:simpleGraph.adjList.keySet())
			{
				for(edgeGroup p:simpleGraph.adjList.get(n))
				{
					if(orientMap.get(n).get(p.otherNode(n))==+1)
						orient.put(p,p.towards(p.otherNode(n)));
					else if(orientMap.get(n).get(p.otherNode(n))==-1)
						orient.put(p,p.towards(n));
					else
						orient.put(p,0);
				}
			}
			orientations.add(orient);
		}
		orientations=checkAndUpdate(orientations);
		return orientations;
	}
	
	/**return all orientations of the simple graph that have incoming flows to all non-sources*/
	private ArrayList<HashMap<edgeGroup,Integer>> checkAndUpdate(ArrayList<HashMap<edgeGroup,Integer>> allOrientations)
	{
		ArrayList<HashMap<edgeGroup,Integer>> corrOrientations=new ArrayList<HashMap<edgeGroup,Integer>>();
		for(HashMap<edgeGroup,Integer> orientation:allOrientations)
		{
			Boolean valid=true;
			for(int n:simpleGraph.adjList.keySet())
			{
				if(componentSources.contains(n) || simpleGraph.adjList.get(n).size()==0)
					continue;
				Boolean in=false;
				for(edgeGroup p:simpleGraph.adjList.get(n))
				{
					if(p.orientToNode(orientation.get(p))==n)
					{
						in=true;
						break;
					}
				}
				if(!in)
				{
					valid=false;
					System.out.println("invalid at node "+n);
					break;
				}
			}
			if(valid)
				corrOrientations.add(orientation);
		}
		System.out.println("Throwing "+(allOrientations.size()-corrOrientations.size())+" orientations, "+corrOrientations.size()+" remain");
		return corrOrientations;
	}
	
	/**map from node to index where it is present*/
	private HashMap<Integer,Integer> revOrient(ArrayList<Integer> orient)
	{
		HashMap<Integer,Integer> rev=new HashMap<Integer,Integer>();
		int l=orient.size();
		for(int i=0;i<l;i++)
		{
			rev.put(orient.get(i), i);
		}
		return rev;
	}
	
	/**do permutations of superset of nodes in the simple graph to determine directions on edges. 
	 * extra nodes are introduced to represent nodes which were excluded in the reduction
	 * return the modified simpleGraph*/
	private HashMap<Integer,Set<Integer>> simpPermute(ArrayList<ArrayList<Integer>> orientations,HashMap<Integer,Integer> extraNode)
	{//list of orientations. each orientation is a list of nodes, in topological sorted order
		HashMap<Integer,Set<Integer>> simpleGraph=toDerivedSimple(extraNode);
		Set<ArrayList<Integer>> orientationSet=new HashSet<ArrayList<Integer>>();
		//nodes will be grouped according to connected components
		//groupLoc gives indices of last element of each group
		
		
		//number of different forests to generate at each step
		for(int i=0;i<maxOrientAtEachStep;i++)
		{
			ArrayList<Integer> or=new ArrayList<Integer>();
			HashMap<Integer,Boolean> visit=new HashMap<Integer,Boolean>();
			HashMap<Integer,ArrayList<Integer>> level=new HashMap<Integer,ArrayList<Integer>>();
			for(int n: simpleGraph.keySet())
			{
				visit.put(n, false);
			}
			
			//only nodes connected to componentSources matter, others dont have any edges connected to them
			for(int n:componentSources)
			{
				if(!visit.get(n))
				{
					level.put(0, new ArrayList<Integer>());
//					level.get(0).add(n);
					simpleDFS(simpleGraph,n,visit,level,0);
					for(int lev=0;;lev++)
					{
						if(!level.containsKey(lev))
							break;
						for(int nodeNo:level.get(lev))
							or.add(nodeNo);
					}
				}
			}
			orientationSet.add(or);
		}
		for(ArrayList<Integer> or:orientationSet)
			orientations.add(or);
//		ArrayList<Integer> groupedNodes=new ArrayList<Integer>(),groupLoc=new ArrayList<Integer>();
//		HashMap<Integer,Boolean> visit=new HashMap<Integer,Boolean>();
//		for(int n: simpleGraph.keySet())
//		{
//			visit.put(n, false);
//		}
//		
//		//only nodes connected to componentSources matter, others dont have any edges connected to them
//		for(int n:componentSources)
//		{
//			if(!visit.get(n))
//			{
//				simpleDFS(simpleGraph,n,groupedNodes,visit);
//				groupLoc.add(groupedNodes.size()-1);
////				nodesCovered+=nodesInGroup;
//			}
//		}
//		ArrayList<ArrayList<Integer>> blankOrientations=new ArrayList<ArrayList<Integer>>();
//		blankOrientations=permuteGroups(groupedNodes,groupLoc,0,blankOrientations);
//		for(ArrayList<Integer> or:blankOrientations)
//			orientations.add(or);
		return simpleGraph;
	}
	
	/**return all permutations where the groups are permuted within themselves and first element of no group is displaced*/
//	private ArrayList<ArrayList<Integer>> permuteGroups(ArrayList<Integer> groupedNodes,ArrayList<Integer> groupLoc,int groupNo,ArrayList<ArrayList<Integer>> orientationsSoFar)
//	{
//		if(groupNo==groupLoc.size())
//			return orientationsSoFar;
//		int grpSize,start;
//		if(groupNo==0)
//		{
//			grpSize=groupLoc.get(groupNo)+1;
//			start=0;
//		}
//		else
//		{
//			grpSize=groupLoc.get(groupNo)-groupLoc.get(groupNo-1);
//			start=groupLoc.get(groupNo-1)+1;
//		}
//		
//		Integer[] group=new Integer[grpSize];
//		for(int i=start;i<=groupLoc.get(groupNo);i++)
//		{
//			group[i-start]=groupedNodes.get(i);
//		}
//		
//		ArrayList<List<Integer>> grpOrientations=new ArrayList<List<Integer>>();
//		
//		permuteEles(group,grpSize,grpSize,grpOrientations);
//		
//		ArrayList<ArrayList<Integer>> orientationsUpdated=new ArrayList<ArrayList<Integer>>();
//		if(orientationsSoFar.size()==0)
//		{
//			
//			for(List<Integer> o:grpOrientations)
//			{
//				ArrayList<Integer> or=new ArrayList<Integer>();
//				for(int oi:o)
//					or.add(oi);
//				orientationsUpdated.add(or);
//			}
//		}
//		for(ArrayList<Integer> oldOrient:orientationsSoFar)
//		{
//			for(List<Integer> addedOrient:grpOrientations)
//			{
//				ArrayList<Integer> updatedOrient=new ArrayList<Integer>();
//				for(Integer i: oldOrient)
//					updatedOrient.add(i);
//				for(Integer i:addedOrient)
//					updatedOrient.add(i);
//				
//				orientationsUpdated.add(updatedOrient);
//			}
//		}
//		
//		return permuteGroups(groupedNodes,groupLoc,groupNo+1,orientationsUpdated);
//	}
//	
//	/**all permutations of elements in group using heap's algorithm*/
//	private void permuteEles(Integer[] group,int size,int n,ArrayList<List<Integer>> grpOrientations)
//	{
////		int sizeLeft;
//		if(size==1)
//		{
//			grpOrientations.add(Arrays.asList(group));
//			return;
//		}
//		if(grpOrientations.size()>1 && end==true)
//			return;
//		for(int i=1;i<size;i++)
//		{
//			permuteEles(group,size-1,n,grpOrientations);
//			
//			if(size%2==1)
//			{
//				Integer temp=group[0];
//				group[0]=group[size-1];
//				group[size-1]=temp;
//			}
//			else
//			{
//				Integer temp=group[i];
//				group[i]=group[size-1];
//				group[size-1]=temp;
//			}
//		}
//	}
	
	/**group nodes by connected components - not doing this now perhaps
	 * do a dfs starting at curr_node, assign each child of curr_node a level of this_level+1
	 * do the recursive call in a random order on each child so that different trees could be formed*/
	private void simpleDFS(HashMap<Integer,Set<Integer>> simpleGraph,int curr_node,HashMap<Integer,Boolean> visit,HashMap<Integer,ArrayList<Integer>> level,int this_level)
	{
		
		visit.put(curr_node, true);
		if(level.containsKey(this_level))
			level.get(this_level).add(curr_node);
		else
		{
			level.put(this_level,new ArrayList<Integer>());
			level.get(this_level).add(curr_node);
		}
		ArrayList<Integer> neighbour=new ArrayList<Integer>();
		for(int n:simpleGraph.get(curr_node))
		{
			if(!visit.get(n))
				neighbour.add(n);
		}
		Collections.shuffle(neighbour);
		for(int n:neighbour)
		{
			simpleDFS(simpleGraph,n,visit,level,this_level+1);
		}
		
		
//		groupedNodes.add(curr_node);
//		visit.put(curr_node, true);
////		noNodesInConnComp+=1;
//		for(int n:simpleGraph.get(curr_node))
//		{
//			if(!visit.get(n))
//			{
//				simpleDFS(simpleGraph,n,groupedNodes,visit);
//			}
//		}
	}
	
	/**if conversion to simple graph and then reducing will give a graph with lesser edges. 
	 * first make the simple graph and bridge set then check for no bridges and no <=2 degree nodes*/
	private Boolean reducible()
	{
//		ArrayList<Integer> reducedSources=new ArrayList<Integer>();
		toSimple();
//		HashMap<Integer,Set<Integer>> simpleAdjList=null;
		bridgeSet();
		Boolean noTwoDeg=true;
		for(int n:simpleGraph.adjList.keySet())
		{
			if(simpleGraph.adjList.get(n).size()<=2 && simpleGraph.adjList.get(n).size()>0 && !componentSources.contains(n))
			{
				noTwoDeg=false;
				break;
			}
		}
		if(noTwoDeg && bridges.size()==0)
			return false;
		else
			return true;
	}
	
	/**dfs for making the reduced graph. 
	 * thisEdge is the edge of the reduced graph, which will either end here, or be continued with the degree is just 2.
	 * To find the degree we remove the bridges, and do the procedure on the separated components.
	 * remEdge is the set of edges that have already been encountered as bridges, so that they are not encountered again during dfs
	 * they are not removed from the adjList to avoid introducing some error. 
	 * returns the id of the highest edge we have so far*/
	private int redGraphDFS(redGraph postReduction,edgeGroup thisEdge,int curr_node,HashMap<Integer,Boolean> visit,Boolean notsource,int highest_edge,HashMap<Integer,Integer> level,Set<edgeGroup> remEdge,Set<edgeGroup> bridges)
	{
		int nextNode;
		visit.put(curr_node, true);
		int degree=0;
		
		//find out the non-bridge adjacency of this
		for(edgeGroup p:simpleGraph.adjList.get(curr_node))
		{
			nextNode=p.otherNode(curr_node);
			
			if(!remEdge.contains(p))
			{
				if(bridges.contains(p))
				{//remove the edge, reduce the next node first
					remEdge.add(p);
					level.put(nextNode, level.get(curr_node)+1);
//					edgeGroup newEdge=new edgeGroup(highest_edge+1,nextNode);
					highest_edge=redGraphDFS(postReduction,null,nextNode,visit,false,highest_edge+1,level,remEdge,bridges);
				}
				else
					degree+=1;
			}
		}
		if(degree==1 && notsource)
		{//no need for this as the edge is a bridge in that case and would have been identified already. this is never encountered
			System.out.println("should not have been encountered!!!!!");
			thisEdge.addEndNode(curr_node);
			postReduction.adjList.get(curr_node).add(thisEdge);
		}
		else if(degree==2 && notsource)
		{
			for(edgeGroup p:simpleGraph.adjList.get(curr_node))
			{
				if(remEdge.contains(p))
					continue;
				
				nextNode=p.otherNode(curr_node);
				
				thisEdge.addEdge(p);
				
				if(!visit.get(nextNode))
				{//tree edge to new vertex
					level.put(nextNode,level.get(curr_node)+1);
					highest_edge=redGraphDFS(postReduction,thisEdge,nextNode,visit,true,highest_edge,level,remEdge,bridges);
				}
				else
				{//this only has 2 edges, so this is either parent or back edge
					if(level.get(curr_node)>level.get(nextNode)+1)
					{//back edge
						if(!postReduction.adjList.get(nextNode).contains(thisEdge))//to prevent adding twice
							postReduction.adjList.get(nextNode).add(thisEdge);
						thisEdge.addEndNode(nextNode);
					}
				}
			}
		}
		else
		{
			for(edgeGroup p:simpleGraph.adjList.get(curr_node))
			{
				if(remEdge.contains(p))
					continue;
				
				nextNode=p.otherNode(curr_node);
				
				if(!visit.get(nextNode))
				{
					highest_edge+=1;
					edgeGroup newEdge=new edgeGroup(highest_edge,curr_node);
					newEdge.addEdge(p);
					postReduction.adjList.get(curr_node).add(newEdge);
					level.put(nextNode,level.get(curr_node)+1);
					highest_edge=redGraphDFS(postReduction,newEdge,nextNode,visit,true,highest_edge,level,remEdge,bridges);
				}
				else
				{
					if(level.get(curr_node)==level.get(nextNode)+1)
					{//parent-child
						postReduction.adjList.get(curr_node).add(thisEdge);
						thisEdge.addEndNode(curr_node);
					}
					else if(level.get(curr_node)>level.get(nextNode)+1)
					{//back edge
//						assert((level.get(curr_node)>level.get(nextNode)+1));
						highest_edge+=1;
						edgeGroup newEdge=new edgeGroup(highest_edge,curr_node,nextNode);
						newEdge.addEdge(p);
						postReduction.adjList.get(curr_node).add(newEdge);
						postReduction.adjList.get(nextNode).add(newEdge);
					}
					else
					{//already taken care of
						
					}
				}
			}
		}
		return highest_edge;
	}
}
