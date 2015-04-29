package io;

import java.io.FileNotFoundException;
import java.io.PrintWriter;

import parameter.Setting;

public class FileOutput {
	private PrintWriter writer;
	private String fileName;
	
	public FileOutput(Setting rc, String label, String head) throws FileNotFoundException{
		fileName = rc.getOutput() + "." + label;
		writer = new PrintWriter(fileName);
		writer.write(head);
		System.out.println("Output: "+fileName);
	}
	
	public FileOutput(Setting rc, String label) throws FileNotFoundException{
		fileName = rc.getOutput() + "." + label;
		writer = new PrintWriter(fileName);
		System.out.println("Output: "+fileName);
	}

	public String getFileName(){
		return fileName;
	}
	
	public void writeline(String line){
		writer.write(line);
	}
	
	public void write(Object item){
		writer.write(item.toString());
		writer.flush();		
	}
	
	public void writemore(Object[] items){
		for(Object item: items){
			writer.write(item.toString());
		}
		writer.flush();		
	}

	public void close(){
		writer.close();
		System.out.println("Close: "+fileName);
	}
}
