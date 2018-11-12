package training;

import com.sun.org.apache.bcel.internal.classfile.Attribute;

import weka.attributeSelection.AttributeSelection;
import weka.attributeSelection.InfoGainAttributeEval;
import weka.attributeSelection.Ranker;
import weka.core.Instances;
import weka.core.Utils;
import weka.core.converters.ConverterUtils.DataSource;
import weka.filters.Filter;
import weka.filters.unsupervised.attribute.NumericToBinary;
import weka.filters.unsupervised.attribute.Remove;

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
		
	}
}
