package training;

import java.util.ArrayList;
import java.util.Enumeration;
import java.util.Random;

import javax.swing.JFrame;

import com.sun.org.apache.bcel.internal.classfile.Attribute;
import com.sun.xml.internal.bind.v2.model.impl.ArrayInfoImpl;

import weka.attributeSelection.AttributeSelection;
import weka.attributeSelection.InfoGainAttributeEval;
import weka.attributeSelection.Ranker;
import weka.classifiers.Classifier;
import weka.classifiers.Evaluation;
import weka.classifiers.evaluation.ThresholdCurve;
import weka.classifiers.trees.J48;
import weka.core.Instances;
import weka.core.Utils;
import weka.core.converters.ConverterUtils.DataSource;
import weka.filters.Filter;
import weka.filters.unsupervised.attribute.NumericToBinary;
import weka.filters.unsupervised.attribute.Remove;
import weka.gui.treevisualizer.PlaceNode2;
import weka.gui.treevisualizer.TreeVisualizer;
import weka.gui.visualize.PlotData2D;
import weka.gui.visualize.ThresholdVisualizePanel;

public class Classification {
	public static void main(String[] args) throws Exception {
		
		/*
		 * Load the data
		 */
		DataSource source = new DataSource("./learningData/training_small.csv");
		Instances data = source.getDataSet();
		System.out.println(data.numInstances() + " instances loaded.");
		

		
		/*
		 * set class attribute numeric to binary
		 */
		NumericToBinary transfer=new NumericToBinary();
		transfer.setOptions(new String[] {"-R",data.numAttributes()+""});
		transfer.setInputFormat(data);
		data=Filter.useFilter(data, transfer);
		
//		System.out.println(data.attributeStats(data.numAttributes()-1));
		data.setClass(data.attribute(data.numAttributes()-1));
//		System.out.println(data.classAttribute());
		
		/*
		 * Feature selection
		 */
		AttributeSelection attSelect = new AttributeSelection();
		InfoGainAttributeEval eval = new InfoGainAttributeEval();
		Ranker search = new Ranker();
		attSelect.setEvaluator(eval);
		attSelect.setSearch(search);
		attSelect.SelectAttributes(data);
		int[] indices = attSelect.selectedAttributes();
		
		System.out.println("Selected attributes: "+Utils.arrayToString(indices));
		for(int index:indices) {
			System.out.print(data.attribute(index).name()+" ");
		}
		
		System.out.println(attSelect.toResultsString());
		
		
		Instances data1=attSelect.reduceDimensionality(data);
		
		
		/*
		 * Build a decision tree
		 */
		String[] options = new String[1];
		options[0] = "-U";
		J48 tree = new J48();
		tree.setOptions(options);
		tree.buildClassifier(data1);
		System.out.println(tree);
		
		
		/*
		 * Visualize decision tree
		 */
		TreeVisualizer tv = new TreeVisualizer(null, tree.graph(),
				new PlaceNode2());
		JFrame frame = new javax.swing.JFrame("Tree Visualizer");
		frame.setSize(800, 500);
		frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		frame.getContentPane().add(tv);
		frame.setVisible(true);
		tv.fitToScreen();
		
		/*
		 * Evaluation
		 */

		Classifier cl = new J48();
		Evaluation eval_roc = new Evaluation(data1);
		eval_roc.crossValidateModel(cl, data1, 10, new Random(1), new Object[] {});
		System.out.println(eval_roc.toSummaryString());
		// Confusion matrix
		double[][] confusionMatrix = eval_roc.confusionMatrix();
		System.out.println(eval_roc.toMatrixString());
		
		/*
		 * Bonus: Plot ROC curve
		 */

		ThresholdCurve tc = new ThresholdCurve();
		int classIndex = 0;
		Instances result = tc.getCurve(eval_roc.predictions(), classIndex);
		// plot curve
		ThresholdVisualizePanel vmc = new ThresholdVisualizePanel();
		vmc.setROCString("(Area under ROC = " + tc.getROCArea(result) + ")");
		vmc.setName(result.relationName());
		PlotData2D tempd = new PlotData2D(result);
		tempd.setPlotName(result.relationName());
		tempd.addInstanceNumberAttribute();
		// specify which points are connected
		boolean[] cp = new boolean[result.numInstances()];
		for (int n = 1; n < cp.length; n++)
			cp[n] = true;
		tempd.setConnectPoints(cp);

		// add plot
		vmc.addPlot(tempd);
		// display curve
		JFrame frameRoc = new javax.swing.JFrame("ROC Curve");
		frameRoc.setSize(800, 500);
		frameRoc.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		frameRoc.getContentPane().add(vmc);
		frameRoc.setVisible(true);
	}
}
