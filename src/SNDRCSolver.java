import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.Paths;
import java.util.*;

import org.jorlib.frameworks.columnGeneration.branchAndPrice.AbstractBranchAndPrice;
import org.jorlib.frameworks.columnGeneration.branchAndPrice.AbstractBranchCreator;
import org.jorlib.frameworks.columnGeneration.branchAndPrice.bapNodeComparators.BFSbapNodeComparator;
import org.jorlib.frameworks.columnGeneration.io.SimpleDebugger;
import org.jorlib.frameworks.columnGeneration.master.cutGeneration.CutHandler;
import org.jorlib.frameworks.columnGeneration.pricing.AbstractPricingProblemSolver;
import org.jorlib.frameworks.columnGeneration.util.Configuration;
import org.jorlib.frameworks.columnGeneration.util.MathProgrammingUtil;

import bap.BranchAndPrice;
import bap.BranchAndPriceA;
import bap.BranchAndPriceA_M;
import bap.BranchAndPriceB;
import bap.BranchAndPriceB_M;
import bap.bapNodeComparators.NodeBoundbapNodeComparator;
import bap.bapNodeComparators.NodeBoundbapNodeComparatorMaxBound;
import bap.branching.BranchOnHoldingEdge;
import bap.branching.BranchOnLocalService;
import bap.branching.BranchOnLocalServiceForAllPricingProblems;
import bap.branching.BranchOnQVarible;
import bap.branching.BranchOnServiceEdge;
import bap.branching.BranchOnServiceEdgeForAllPricingProblems;
import bap.branching.BranchOnTimeService;
import bap.branching.BranchOnTimeServiceForAllPricingProblems;
import cg.Cycle;
import cg.ExactPricingProblemSolver;
import cg.SNDRCPricingProblem;
import cg.master.Master;
import cg.master.SNDRCMasterData;
import cg.master.cuts.StrongInequalityGenerator;
import ch.qos.logback.core.net.SyslogOutputStream;
import logger.BapLogger;
import logger.BapLoggerA;
import logger.BapLoggerA_M;
import logger.BapLoggerB;
import logger.BapLoggerB_M;
import model.SNDRC;
import model.SNDRC.Edge;
import sun.security.krb5.PrincipalName;

public class SNDRCSolver {
    SNDRC dataModel;
    boolean ifOptGetFromSubGraph;
    BranchAndPriceA bap;

    public SNDRCSolver(SNDRC dataModel,String fileName) throws IOException {
        this.dataModel = dataModel;

        // Create the pricing problems
        List<SNDRCPricingProblem> pricingProblems = new LinkedList<SNDRCPricingProblem>();
        for (int capacityType = 0; capacityType < dataModel.numOfCapacity; capacityType++) {
            for (int originNode = 0; originNode < dataModel.numNode; originNode++) {
                String name = "capacity type: " + capacityType + " origin node: " + originNode;
                SNDRCPricingProblem pricingProblem = new SNDRCPricingProblem(dataModel, name, capacityType, originNode);
                pricingProblems.add(pricingProblem);
            }
        }

        // Create a cutHandler
        CutHandler<SNDRC, SNDRCMasterData> cutHandler = new CutHandler<>();
        StrongInequalityGenerator cutGen = new StrongInequalityGenerator(dataModel, pricingProblems, 0);
        // cutHandler.addCutGenerator(cutGen);

        // Create the Master Problem
        Master master = new Master(dataModel, pricingProblems, cutHandler, cutGen, false);

        // Define which solvers to use
        List<Class<? extends AbstractPricingProblemSolver<SNDRC, Cycle, SNDRCPricingProblem>>> solvers = Collections
                .singletonList(ExactPricingProblemSolver.class);

        // Define one or more Branch creators
        // List<? extends AbstractBranchCreator<SNDRC, Cycle,
        // SNDRCPricingProblem>> branchCreators=Arrays.asList(new
        // BranchOnQVarible(dataModel, pricingProblems,0.5),new
        // BranchOnServiceEdgeForAllPricingProblems(dataModel, pricingProblems,
        // 0.5),new BranchOnServiceEdge(dataModel, pricingProblems, 0.5));
        List<? extends AbstractBranchCreator<SNDRC, Cycle, SNDRCPricingProblem>> branchCreators = Arrays.asList(
                new BranchOnLocalServiceForAllPricingProblems(dataModel, pricingProblems, 0.5),
                new BranchOnLocalService(dataModel, pricingProblems, 0.5),
                new BranchOnServiceEdge(dataModel, pricingProblems, 0.5));
        // List<? extends AbstractBranchCreator<SNDRC, Cycle,
        // SNDRCPricingProblem>> branchCreators=Arrays.asList(new
        // BranchOnTimeServiceForAllPricingProblems(dataModel, pricingProblems,
        // 0.5, 5),new BranchOnServiceEdge(dataModel, pricingProblems, 0.5));

        // Create a Branch-and-Price instance
         this.bap=new BranchAndPriceA(dataModel, master,pricingProblems, solvers,branchCreators,Double.MAX_VALUE,0.6,0.2,0.1,10,0.0001,3,false);
        // BranchAndPriceB bap=new BranchAndPriceB(dataModel, master,
        // pricingProblems, solvers,
        // branchCreators,Double.MAX_VALUE,0.65,0.2,0.1,1,0.001,3,0.1,true);
//        BranchAndPriceB_M bap = new BranchAndPriceB_M(dataModel, master, pricingProblems, solvers, branchCreators,Double.MAX_VALUE, 0.65, 0.3, 0.1, 1, -0.001,3, 0, true, false);
//        BranchAndPriceB_M bap = new BranchAndPriceB_M(dataModel, master, pricingProblems, solvers, branchCreators,5000000, 0.65, 0.3, 0.1, 1, -0.001,3, 0, true, false);
        // BranchAndPriceA_M bap=new BranchAndPriceA_M(dataModel, master,
        // pricingProblems, solvers,
        // branchCreators,Double.MAX_VALUE,0.6,0.3,0.1,10,0.001,10,0.1,false,true);
//        BranchAndPriceA_M bap = new BranchAndPriceA_M(dataModel, master, pricingProblems, solvers, branchCreators,5000000, 0.6, 0.2, 0.5, 10, 0.001, 10, 0.1, false, true);
        // bap.setNodeOrdering(new BFSbapNodeComparator());
        // bap.setNodeOrdering(new NodeBoundbapNodeComparatorMaxBound());
        bap.setNodeOrdering(new NodeBoundbapNodeComparator());

        // OPTIONAL: Attach a debugger
         SimpleDebugger debugger=new SimpleDebugger(bap, true);

        // OPTIONAL: Attach a logger to the Branch-and-Price procedure.
        // BapLoggerA logger=new BapLoggerA(bap, new
        // File("./output/BAPlogger.log"));
        // BapLoggerB logger=new BapLoggerB(bap, new
        // File("./output/BAPlogger.log"));
//        BapLoggerB_M logger = new BapLoggerB_M(bap, new File("./output/BAPlogger.log"));
//        BapLoggerB_M logger = new BapLoggerB_M(bap, new File("./output/"+fileName+".log"));
//        BapLoggerA_M logger = new BapLoggerA_M(bap, new File("./output/BAPlogger.log"));
        BapLoggerA logger = new BapLoggerA(bap, new File("./output/BAPlogger.log"));

        // Solve the TSP problem through Branch-and-Price
        // bap.runBranchAndPrice(System.currentTimeMillis()+18000000L); //5
        // hours
        bap.runBranchAndPrice(System.currentTimeMillis() + 7200000L); // 2
                                                                       // hours

        this.ifOptGetFromSubGraph = bap.GetIfOptGetFromSubGraph();

        // Print solution:
        System.out.println("================ Solution ================");
        System.out.println("BAP terminated with objective : " + bap.getObjective());
        System.out.println("Total Number of iterations: " + bap.getTotalNrIterations());
        System.out.println("Total Number of processed nodes: " + bap.getNumberOfProcessedNodes());
        System.out.println("Total Time spent on master problems: " + bap.getMasterSolveTime()
                + " Total time spent on pricing problems: " + bap.getPricingSolveTime());
        System.out.println("Best bound : " + bap.getBound());

        // Count the number of key edges
        Set<Integer> keyServiceEdgeIndexSet = new TreeSet<>();
        if (bap.hasSolution()) {
            List<Cycle> solution = bap.getSolution();
            for (Cycle column : solution) {
                Set<Integer> tempEdgeSet = column.edgeIndexSet;
                for (int edgeIndex : tempEdgeSet) {
                    if (((ifOptGetFromSubGraph) && (dataModel.subEdgeSet.get(edgeIndex).edgeType == 0))
                            || ((!ifOptGetFromSubGraph) && (dataModel.edgeSet.get(edgeIndex).edgeType == 0))) {
                        
                        int index;
                        if(ifOptGetFromSubGraph){
                            index=dataModel.edgeSetIndexMap.get(edgeIndex);
                        }else{
                            index=edgeIndex;
                        }
                        
                        if (!keyServiceEdgeIndexSet.contains(index)) {
                            keyServiceEdgeIndexSet.add(index);
                        }
                    }
                }
            }
        }

        System.out.println(keyServiceEdgeIndexSet.toString());
        System.out.println("The number of service edges used= " + keyServiceEdgeIndexSet.size());
        System.out.println();

        // record the information of costs on vehicle and commodity and calculate "no load ratio"
        int vehicleFixTotalCost = 0;
        Map<Cycle, Integer> vehicleVarCost = new HashMap<>();
        int vehicleVarTotalCost = 0;
        List<Double> commodityCost = new ArrayList<>();
        Double commodityTotalCost = (double) 0;
        
        long vehicleDowork=0;
        double commodityDowork=0;
        
        int totalNumVehicle=0;
        Map<Integer,Integer> vehicleCoverServiceEdgeRecord=new HashMap<>();
        double[][] commodityFlowIntoTerminal=new double[dataModel.numNode][dataModel.timePeriod];

        if (bap.hasSolution()) {
            System.out.println("Solution is optimal: " + bap.isOptimal());
            System.out.println("Columns (only non-zero columns are returned):");
            List<Cycle> solution = bap.getSolution();
            for (Cycle column : solution) {
                System.out.println(column);
                System.out.println(out(column) + ":" + bap.GetOptSolutionValueMap().get(column));
                totalNumVehicle+=MathProgrammingUtil.doubleToInt((double) bap.GetOptSolutionValueMap().get(column));

                int fixCost = (int) dataModel.fixedCost[column.associatedPricingProblem.originNodeO][column.associatedPricingProblem.capacityTypeS];
                int varCost = (int) (column.cost - fixCost);

                double columnValue= (double) bap.GetOptSolutionValueMap().get(column);
                int value=MathProgrammingUtil.doubleToInt(columnValue);
//                column.value=value;
                
                vehicleFixTotalCost += fixCost*columnValue;
                vehicleVarTotalCost += varCost*columnValue;
                vehicleVarCost.put(column, varCost);
                System.out.println("Fix cost= " + fixCost + " variable cost= " + varCost);

                System.out.println();
                
                
                for(int edgeIndex:column.edgeIndexSet){
                    
                	int index=edgeIndex;
                    Edge edge;
                    if(!ifOptGetFromSubGraph){
                        edge=dataModel.edgeSet.get(edgeIndex);
                    }else {
                    	edge=dataModel.subEdgeSet.get(edgeIndex);
                    	index=dataModel.edgeSetIndexMap.get(edgeIndex);
                    }
                    if(edge.edgeType==0){
                        vehicleDowork+=dataModel.capacity[column.associatedPricingProblem.capacityTypeS]*edge.duration*value; 
                        
                        if(!vehicleCoverServiceEdgeRecord.containsKey(index)) {
                        	vehicleCoverServiceEdgeRecord.put(index, value);
                        }else {
                        	vehicleCoverServiceEdgeRecord.put(index, vehicleCoverServiceEdgeRecord.get(index)+value);
                        }
                    }
                    
                    
                    
                }
            }

        }

        // calculate cost of commodities
        List<Map<Integer, Double>> optXValues = bap.GetOptXValues();
        for (int commodity = 0; commodity < dataModel.numDemand; commodity++) {
            Map<Integer, Double> xValues = optXValues.get(commodity);
            double cost = 0;
            for (int edgeIndex : xValues.keySet()) {
                Edge edge;
                if (!ifOptGetFromSubGraph) {
                    edge = dataModel.edgeSet.get(edgeIndex);
                } else {
                    edge = dataModel.subEdgeSet.get(edgeIndex);
                }

                if (edge.edgeType == 0) {
                    cost += dataModel.beta * edge.duration * xValues.get(edgeIndex);
                    
                    commodityDowork+=edge.duration*xValues.get(edgeIndex);
                    
//                    commodityFlowIntoTerminal[edge.v][edge.t2]+=MathProgrammingUtil.doubleToInt(xValues.get(edgeIndex));
                    commodityFlowIntoTerminal[edge.v][edge.t2]+=xValues.get(edgeIndex);
                }
            }

            commodityCost.add(cost);
            commodityTotalCost += cost;
        }

        System.out.println("fix cost+variable cost+commodity cost= " + vehicleFixTotalCost + "+" + vehicleVarTotalCost
                + "+" + commodityTotalCost + "=" + (vehicleFixTotalCost + vehicleVarTotalCost + commodityTotalCost));
        System.out.println();
        System.out.println("vehicle dowork= "+vehicleDowork+" commodity dowork= "+commodityDowork);
        System.out.println("no load ratio= "+(vehicleDowork-commodityDowork)/vehicleDowork);
        System.out.println();
        System.out.println("Total vehicles used= "+totalNumVehicle);
        System.out.println();
        System.out.println("vehicleCoverServiceEdge information:");
        System.out.println(vehicleCoverServiceEdgeRecord.toString());
        System.out.println();
        System.out.println("commodityFlowIntoTerminal information:");
        for(int terminal=0;terminal<dataModel.numNode;terminal++) {
        	for(int i=0;i<dataModel.timePeriod;i++) {
        		System.out.print((int)commodityFlowIntoTerminal[terminal][i]+" ");
        	}
        	System.out.println();
//        	System.out.println(Arrays.toString(commodityFlowIntoTerminal[terminal]));
        }
        
        
        System.out.println();
        System.out.println("vehicle pattern information:");
        Set<int[]> avoidRepeatPattern=new HashSet<>();
        List<Cycle> solution = bap.getSolution();
        for (Cycle column : solution) {
        	boolean check=true;
        	for(int[] recordPattern:avoidRepeatPattern) {
        		if(Arrays.equals(recordPattern, column.pattern)) {
        			check=false;
        			break;
        		}
        	}
        	
        	if(check) {
        		System.out.println(Arrays.toString(column.pattern));
        		avoidRepeatPattern.add(column.pattern);
        	}
        }
        
        
//        
//        //compare keyServiceEdgeIndexSet and serviceEdgeSet0
//        Scanner in = new Scanner(Paths.get("./data/testset/test12D_compareInfo.txt"));
//        String line=in.nextLine();
//        String[] result=line.split(", ");
//        Set<Integer> serviceEdgeSet0=new HashSet<>();
//        for(int i=0;i<result.length;i++){
//            serviceEdgeSet0.add(Integer.parseInt(result[i]));
//        }
//        
//        int commonCount=0;
//        for(int edgeIndex:keyServiceEdgeIndexSet){
//            if(serviceEdgeSet0.contains(edgeIndex)){
//                commonCount++;
//            }
//        }
//        
//        System.out.println("Common edge information:");
//        System.out.println("size0= "+serviceEdgeSet0.size()+" size1= "+keyServiceEdgeIndexSet.size()+" # of common edges= "+commonCount);
//        double ratio=(double)commonCount/serviceEdgeSet0.size();
//        System.out.println("common edge ratio= "+ratio);
//        
//        
//        //compare avoidRepeatPattern and pattern0
//        Set<int[]> pattern0=new HashSet<>();
//        while(in.hasNextLine()){
//            line=in.nextLine();
//            result=line.split(", ");
//            int[] temp=new int[result.length];
//            for(int i=0;i<result.length;i++){
//                temp[i]=Integer.parseInt(result[i]);
//            }
//            pattern0.add(temp);
//        }
//        
//        commonCount=0;
//        for(int[] pattern:avoidRepeatPattern){
//            boolean check=false;
//            for(int[] ele:pattern0){
//                if(Arrays.equals(ele, pattern)){
//                    check=true;
//                    break;
//                }
//            }
//            
//            if(check) commonCount++;
//        }
//        
//        System.out.println("Common pattern information:");
//        System.out.println("size of pattrenSet0= "+pattern0.size()+" size of patternSet1= "+avoidRepeatPattern.size()+" # of common patterns= "+commonCount);
//        ratio=(double)commonCount/pattern0.size();
//        System.out.println("common pattern ratio= "+ratio);
        
        
        
        

        // output x variables
        for (int demand = 0; demand < dataModel.numDemand; demand++) {
            for (int edgeIndex : optXValues.get(demand).keySet()) {
                if (optXValues.get(demand).get(edgeIndex) > 0.01) {
                    Edge edge;

                    if (!ifOptGetFromSubGraph) {
                        edge = dataModel.edgeSet.get(edgeIndex);
                    } else {
                        edge = dataModel.subEdgeSet.get(edgeIndex);
                    }

                    
                    if(edge.edgeType==0){
                        System.out.println("x[" + demand + "]:" + edge.u + "," + edge.t1 + "->" + edge.v + "," + edge.t2
                                + "= " + optXValues.get(demand).get(edgeIndex) + " " + edge.duration);
                    }

                }
            }
            System.out.println("total cost= " + commodityCost.get(demand));
            System.out.println();
        }

        bap.close();
        cutHandler.close();

    }

    public String out(Cycle column) {

        Queue<Edge> path = new PriorityQueue<>();

        if (!ifOptGetFromSubGraph) {
            for (int edgeIndex : column.edgeIndexSet) {
                path.add(dataModel.edgeSet.get(edgeIndex));
            }
        } else {
            for (int edgeIndex : column.edgeIndexSet) {
                path.add(dataModel.subEdgeSet.get(edgeIndex));
            }
        }

        StringBuilder pathRecord = new StringBuilder();

        // int count=0;
        // for(Edge edge:path) {
        // count++;
        // pathRecord.append(edge.start);
        //
        // if(count!=column.edgeIndexSet.size()) {
        // pathRecord.append("->");
        // }
        //
        // }

        Edge edge = null;
        int size = path.size();
        for (int i = 0; i < size; i++) {

            edge = path.poll();
            pathRecord.append('(');
            pathRecord.append(edge.u);
            pathRecord.append(',');
            pathRecord.append(edge.t1);
            pathRecord.append(')');

            pathRecord.append("->");

        }

        pathRecord.append('(');
        pathRecord.append(edge.v);
        pathRecord.append(',');
        pathRecord.append(edge.t2);
        pathRecord.append(')');

        return pathRecord.toString();

    }

    
    public void output(String filename,BranchAndPriceA bap) throws FileNotFoundException{
    	PrintWriter out=new PrintWriter(filename);
    	
    	out.println("LP information:");
    	for(int k=0;k<dataModel.numDemand;k++){
    		Map<Integer,Double> map=(Map<Integer, Double>) bap.GetxValuesForRootLP().get(k);
    		out.println(k);
    		for(int edgeIndex:map.keySet()){
    			double value=map.get(edgeIndex);
    			if(value>0.0001&&dataModel.edgeSet.get(edgeIndex).edgeType==0){
    				out.print(edgeIndex+" "+value+" ");
    			}
    		}
    		out.println();
    		
    	}
    	
    	//Here we need to notice the change of edge index
    	out.println();
    	out.println("Final Solution:");
    	for(int k=0;k<dataModel.numDemand;k++){
    		Map<Integer,Double> map=(Map<Integer, Double>) bap.GetOptXValues().get(k);
    		out.println(k);
    		
    		for(int edgeIndex:map.keySet()){
    			double value=map.get(edgeIndex);
    			
    			int index;
    			if(ifOptGetFromSubGraph){
    				index=dataModel.edgeSetIndexMap.get(edgeIndex);
    			}else index=edgeIndex;
    			
    			if(value>0.0001&&dataModel.edgeSet.get(index).edgeType==0){
                    out.print(index+" "+value+" ");
    			}
    			
    		}
    		out.println();
    		
    	}
    	
    	out.close();
    	
    }
    
	public static void main(String[] args) throws IOException {
    	
    	SNDRC sndrc=new SNDRC("./data/testset/example.txt");
//    	sndrc.outputFeature("./learningData/result/1-1.txt");
    	SNDRCSolver solver=new SNDRCSolver(sndrc,"BAPlogger");
//    	solver.output("./learningData/result/1-2.txt", solver.bap);
    	
//    	int count=0;
//    	String path="../learningData/test";
//    	File file=new File(path);
//    	File[] fs=file.listFiles();
//    	
//    	for(File f:fs){
//    		if(!f.isDirectory()&&!f.isHidden()){
//    			String filename=count+"_0.txt";
//    			
//    		    System.out.println("Sovle for "+f.toString());
//    		    System.out.println();
//    	    	SNDRC sndrc=new SNDRC(f.toString());
//    	    	sndrc.outputFeature("../learningData/result/"+filename);
//    	    	SNDRCSolver solver=new SNDRCSolver(sndrc,"BAPlogger");
//    	    	filename=count+"_1.txt";
//    	    	solver.output("../learningData/result/"+filename, solver.bap);
//    	    	
//    	    	count++;
//    			System.out.println();
//    		}
//    	}

//        SNDRC sndrc;
//        for (String arg : args) {
//            long time0 = System.currentTimeMillis();
//            sndrc = new SNDRC(arg);
//            // sndrc.Output();
//            Properties properties = new Properties();
//            // properties.setProperty("EXPORT_MODEL", "True");
//            // properties.setProperty("MAXTHREADS", "10");
//            // properties.setProperty("PRECISION", "0.001");
//            properties.setProperty("CUTSENABLED", "false");
//            Configuration.readFromFile(properties);
//
//            new SNDRCSolver(sndrc,"BAPlogger");
//
//            long time1 = System.currentTimeMillis();
//            System.out.println();
//            System.out.println("Total time= " + (time1 - time0));
//
//        }
        
        
//        SNDRC sndrc;
//        String path="./data/transferData/transfer1/test12D/";
//        String name0="test12D";
//        double[] var={0.5,1.0,2.0,5.0,10.0};
//        
//        Properties properties = new Properties();
//        properties.setProperty("CUTSENABLED", "false");
//        Configuration.readFromFile(properties);
//        
//        int count=0;
//        for(int i=0;i<var.length;i++){
//            double variance=var[i];
//            
//            for(int j=0;j<5;j++){
//                count++;
//                
//                if(count>24){
//                    String name=path+name0+"_"+variance+"_"+j+".txt";
//                    System.out.println("Solve for "+name);
//                    System.out.println();
//                    
//                    
//                    long time0 = System.currentTimeMillis();
//                    sndrc = new SNDRC(name);
//                    
//
//                    
//                    new SNDRCSolver(sndrc,name0+"_"+variance+"_"+j);
//                    
//                    long time1 = System.currentTimeMillis();
//                    System.out.println();
//                    System.out.println("Total time= " + (time1 - time0));
//                    System.out.println();
//                }
//                
//
//                
//            }
//        }
        
        
//        path="./data/transferData/transfer2/test12D/";
//        name0="test12D";
////        double[] times={1.2,1.5,1.8,2.0,3.0};
//        double[] times={1.2};
//        
//        for(double time:times){
//            String name=path+name0+"_"+time+".txt";
//            System.out.println("Solve for "+name);
//            System.out.println();
//            
//            
//            long time0 = System.currentTimeMillis();
//            sndrc = new SNDRC(name);
//            
//
//            
//            new SNDRCSolver(sndrc,name0+"_"+time);
//            
//            long time1 = System.currentTimeMillis();
//            System.out.println();
//            System.out.println("Total time= " + (time1 - time0));
//            System.out.println();
//        }
    	
    	
        
        

              
        

    }

}
