package optimizer;

import java.util.HashMap;
import java.util.HashSet;
import java.util.Set;

public class edgeGroup {
	private final int ID;
	private int node1,node2;
	private Set<edgeGroup> edgeMapping;
	edgeGroup(int id,int node1)
	{
		this.ID=id;
		this.node1=node1;
		this.node2=-1;
		edgeMapping=new HashSet<edgeGroup>();
	}
	edgeGroup(int id,int node1,int node2)
	{
		this.ID=id;
		this.node1=node1;
		this.node2=node2;
		edgeMapping=new HashSet<edgeGroup>();
	}
	
	/**if node1 is given return node2 and vice-versa*/
	int otherNode(int nodei)
	{
		if(nodei==node1)
			return node2;
		else if(nodei==node2)
			return node1;
		else
			return -1;
	}
	
	void addEdge(edgeGroup prevEdge)
	{
		edgeMapping.add(prevEdge);
	}
	
	void addEndNode(int nodei)
	{
		if(this.node1==-1)
			this.node1=nodei;
		else
			this.node2=nodei;
	}
	
	int getID()
	{
		return ID;
	}
	int getMappingSize()
	{
		return edgeMapping.size();
	}
	Set<edgeGroup> getMapping()
	{
		return edgeMapping;
	}
	int getNode1()
	{
		return node1;
	}
	int getNode2()
	{
		return node2;
	}
	int towards(int nodei)
	{
		if(nodei==node1)
			return -1;
		else
		{
			assert(nodei==node2);
			return 1;
		}
			
	}
	int orientToNode(int orientation)
	{
		if(orientation==1)
			return node2;
		else if(orientation==-1)
			return node1;
		else
			return -1;
	}
	
	

}
