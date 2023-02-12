package a_linesOfData;

import java.awt.Color;
import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Vector;

import javax.vecmath.Vector3d;

import java.io.FileReader;
import java.io.BufferedOutputStream;
import java.io.BufferedReader;
import java.io.IOException;
import java.io.FileOutputStream;
import java.io.ObjectOutputStream;

import processing.core.PGraphics3D;
import bsim.BSim;
import bsim.BSimChemicalField;
import bsim.BSimTicker;
import bsim.BSimUtils;
import bsim.export.*;
import bsim.draw.BSimP3DDrawer;
import bsim.export.BSimLogger;
import bsim.particle.BSimBacterium;
import java.util.Random;
/**
 * Tests out the functionality of the {@link BSimChemicalField}.</br>
 * A number of {@link BSimBacterium} bacteria are set up to swim around in a chemical field.
 * The bacteria deposit chemical randomly in space.
 * This acts as an attractor for other bacteria in the vicinity (chemotaxis).
 */


public class lineOfData {

	/*********************************************************
	 * Simulation Definition
	 *********************************************************/
	public static void main(String[] args) {

		/*********************************************************
		 * Create a new simulation object and set up simulation settings
		 */
		

		final int bounds = 128;
		BSim sim = new BSim();		
		sim.setDt(1);//default is 0.01
		sim.setTimeFormat("0.00");
		sim.setSolid(true,true,true);
		sim.setBound(bounds,bounds,bounds);
		sim.setSimulationTime(150);		// Total simulation time [sec]
		
		
		
			
		/*********************************************************
		 * Set up the chemical field
		 */
		boolean average_over_sphere = false;
		final int boxno = bounds;
		final int xyboxno = bounds; //changing boxes sizes and nos changes the graph. can be due to accuracy of the results. change is slight
		final double c = 5e9; // molecules (controls intensity with no major graph difference in morphology)
		final double decayRate = 1e-5; // 4.44e-4 good Can regulate decay rate with pulse to get the top hat shape. 
		final double diffusivity = 0.1; // (microns)^2/sec
		final double threshold = 3e8; //Threshold of FDF chemical 7e8 for cfgl
		final double quantityadd = 5e9; //Quanitity of chemical added
		final double firstwave = 1000; // Time slot where first wave stops
		final double secondwave = 1440; // Time slot for second wave beginning
		final BSimChemicalField field = new BSimChemicalField(sim, new int[]
				{xyboxno,xyboxno,boxno}, diffusivity, decayRate);
//		field.linearZ(0,c);
		
		
		
		
		//Set temporary number and total number of bacteria
		
		int totbac = 512; // Bacterial comes from 33x33 grid
		
		//For even space
		double tempnum = 1/(2*totbac);
		//For packed
		//double tempnum = 1;
		
		/*********************************************************
		 * Set up the bacteria
		 */
		
		BufferedReader reader;
		String positions = null;
		String line = null;
		Vector<String> coordinates = new Vector<String>();
		
		
		// Read in the data file. Shoving algorithm doesn't work so need to make sure points
				// aren't overlapping
				try 
				{
					reader = new BufferedReader(new FileReader("C:\\Users\\user\\OneDrive\\Documents\\Uni\\4th year\\MPhys\\rough shapes\\mengerShapes\\flatMenger.csv"));
					// Read the first line in the file
					line = reader.readLine();

					while (line != null) 
					{
						// Store coordinates in a vector
						coordinates.add(line);
						// get new line
						line = reader.readLine();
						
					}
					
					reader.close();
					System.out.println("File read");
				}
				
				catch(IOException e)
				{
					System.out.println("Problem reading file.");
				}
		
		
		final Vector<BSimBacterium> bacteria = new Vector<BSimBacterium>();		
		while(bacteria.size() < totbac) 
		{
			
			int bacteriaNumber = bacteria.size();
		
		
			////////////////////////////////////////////////////////////////////////////////
			//////////////////// Generate position of bacteria /////////////////////////////
			////////////////////////////////////////////////////////////////////////////////
			
			// get the positions as a string
			positions = coordinates.get(bacteriaNumber);
			
			// Separate the coordinates
			String[] bacPosition = positions.split(",", 3);
			
			String x0 = bacPosition[0];
			String y0 = bacPosition[1];
			String z0 = bacPosition[2];
			
			// Remove all non-integer characters
			x0.replaceAll("[^0-9]+", "");
			y0.replaceAll("[^0-9]+", "");
			z0.replaceAll("[^0-9]+", "");
			
			// Turn string to double.
			double x00 = Double.parseDouble(x0);
			double z00 =  Double.parseDouble(y0); // Swap x and y so shape is right way up
			double y00 = Double.parseDouble(z0);

			double rad = 0.05;
			////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////
			 
			
			BSimBacterium p = new BSimBacterium(sim, 
				//This uses random bacteria	
					
					//Spherical Test
						 	
					new Vector3d(37 +(x00), 
							0.5 * bounds + y00, 
							37 + z00)) {
				
				// This uses the even distancing
					//new Vector3d(/*0*Math.random()*/0.5*sim.getBound().x, 
					//	/*0*Math.random()*/0.5*sim.getBound().y, 
					//	/*0*Math.random()*/0.5*(tempnum)*sim.getBound().z)) {
				// Bacteria move etc. and also add chemical to the global field.
				public void action() {
					super.action();
					/*if (Math.random() < sim.getDt())
						field.addQuantity(position, 1e9);*/
				
					//Condition of Pulse
					
					//if ((sim.getTime() > 60 && sim.getTime() < 64 && (x00/sim.getBound().x < rad && y00/sim.getBound().y < rad && z00/sim.getBound().z < rad && x00/sim.getBound().x > -rad && y00/sim.getBound().y > -rad && z00/sim.getBound().z > -rad))  || (field.getConc(position) > 1e4 && sim.getTime() < 300 && field.getConc(position)<5e9))
						//field.addQuantity(position, 5e9);
						//System.out.println("Fire");
					
					if ( sim.getTime()%1440 > 0 && sim.getTime()%1440 < 4 ) 
					{
							field.addQuantity(64, 64, 37, quantityadd);
						/*
						for (int i =0; i <54; i++)
						{
							field.addQuantity(64, 64, 37 + i, quantityadd);
						}*/

					}
					if ((field.getConc(position) > threshold && 
							  sim.getTime() < firstwave))
							field.addQuantity(position, quantityadd);
						//if (field.getConc(position) < threshold/2 && sim.getTime() > firstwave)
							//field.addQuantity(position, 0.001*quantityadd);
						
						if (field.getConc(position) > threshold && 
							sim.getTime()>secondwave && 
							field.getConc(position) < quantityadd)
							field.addQuantity(position, quantityadd);
						
						
						
					
		//position is the box position and q is the quantity of chemical released in the box. makes curves more prominent up to a certain point		
								
					
					
					/*if (sim.getTime() > 4800 && (sim.getTime()) % 10 == 0 && sim.getTime() < 8400 && Math.random() > 0.99)
						field.addQuantity(position, 5e9)
					;*/
				} 
			};
			// Chemotaxis according to chemical field strength
			p.setGoal(field);
			
			if(!p.intersection(bacteria)) bacteria.add(p);
			
			//Update the temporary number for different position for new bacteria
			//Evenly Spread
			//tempnum = tempnum + 1/totbac;
			
			//Packed
			tempnum = tempnum + 0.015;
		}
		
		
		
		/*********************************************************
		 * Set the ticker, define what happens each time step
		 */
		sim.setTicker(new BSimTicker() {
		
			public void tick() {
				for(BSimBacterium b : bacteria) {
					b.action();		
					//b.updatePosition();
					
				}
				field.update();
				 
			}
			
		});

						
		/*********************************************************
		 * Set the drawer for the simulation
		 */
		sim.setDrawer(new BSimP3DDrawer(sim, 600,600) {
			@Override
			public void scene(PGraphics3D p3d) {	
				draw(field, Color.CYAN, (float)(0.1*255/threshold));						
				for(BSimBacterium p : bacteria) draw(p, Color.GREEN);	  // Colour the bacteria	
			}
		});	
		
		String resultsDir = BSimUtils.generateDirectoryPath("./results/rough_surface/Menger/");
		
		// Name the csv file
			BSimLogger logger = new BSimLogger(sim, resultsDir + "flatMenger.csv") {
			
			@Override
			
			public void before() {
				super.before();
				// Write a header containing the names of variables we will be exporting
				write("time,Intensity"); 
			}
			
		
			public void during() {
                // Counter for the number of current collisions
                //double Intensity = field.getConc(0, 0, 0);
                //List<Double> intensities = new ArrayList<>();
                // Loop through the bacteria and count collisions
                //for (BSimTutorialBacterium p : bacteria){
                //String intensities = "";
                String intensities = "";
                String totchem = "";
                //At infinity the diffusion should be closer to each next distance


                
                if(average_over_sphere){

                    for(int r = 0; r<boxno/2; r++){
                        if(r==0){
                            intensities += "," + (field.getConc(xyboxno/2,xyboxno/2,boxno/2));
                        }else{
                            Vector<Double> intensities_on_sphere = new Vector<Double>();
                            // generate points on a sphere of radius r
                            final int points_density = 10;
                            // determine the number of points on the sphere from the density
                            int points_nb = (int) (points_density* Math.pow(r,2));
                            // generate the random points
                            for(int i=0; i<=points_nb;i++){
                                // x,y,z for points on a sphere
                                double x = new Random().nextGaussian(0,1);
                                double y = new Random().nextGaussian(0,1);
                                double z = new Random().nextGaussian(0,1);
                                // normalise to a sphere
                                double normalisation = r/Math.sqrt(x*x + y*y + z*z);

                                // change units to be box index
                                normalisation *= (float)(boxno)/(float)(bounds);

                                // generate the final x, y and z
                                x *= normalisation;
                                y *= normalisation;
                                z *= normalisation;



                                // for details on these last few steps, refer to https://www.notion.so/Lab-book-01dac4b0eedd4f47b7bee779863b8507#cf3faa9db8fc4f34a6030758018426e8

                                // round the coordinates up and displaces them to the center of the sim, read the intensity at this point on the sphere and store it
                                intensities_on_sphere.add(field.getConc((int)(Math.ceil(x) + boxno/2),
                                                                        (int)(Math.ceil(y) + boxno/2),
                                                                        (int)(Math.ceil(z) + boxno/2)));

                            }

                            // average over all the points in intensities_on_sphere
                            double sum = 0;
                            for(Double v: intensities_on_sphere){
                                sum += v;
                            }
                            // write to file the intensity
                            intensities += "," + sum/ intensities_on_sphere.size();
                            totchem += "," + field.totalQuantity();
                        }

                    }



                }else{ // give the intensities along the z-axis
                    for (int i = 0; i < 54; i++) {
                        //Intensity = field.getConc(0, 0, i);
                    	
                        totchem += "," + field.totalQuantity();
                        
                        
                        intensities += "," + (field.getConc(64, 64, 37 + i)); // these are supposed to be the indexes of the boxes coordinates, ie boxsize dependent
                        
                    }

                }

                write(sim.getFormattedTime()  + intensities);



                //}

                // Write the time and number of collisions to file

            }
        };
        sim.addExporter(logger);
        

        // Run the simulation

        sim.export();
        //sim.preview();



    }

}
