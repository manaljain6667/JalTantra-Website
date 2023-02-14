package optimizer;

import org.apache.commons.math3.analysis.ParametricUnivariateFunction;


public class CostFunc implements ParametricUnivariateFunction {
	public double value(double x, double...parameters )
	{
		return parameters[0]*Math.pow(x,parameters[1])+parameters[2];
	}
	public double[] gradient(double x, double...parameters )
	{
		double g[]=new double[3];
		g[0]=Math.pow(x,parameters[1]);
		g[1]=parameters[0]*Math.pow(x,parameters[1])*Math.log(x);
		g[2]=1;
		
		return g;
	}

}
