package a_planeOfData;

import java.awt.Color;

import java.util.ArrayList;
import java.util.List;
import java.util.Vector;

import javax.vecmath.Vector3d;

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


public class a_planeOfData {

	/*********************************************************
	 * Simulation Definition
	 *********************************************************/
	public static void main(String[] args) {

		/*********************************************************
		 * Create a new simulation object and set up simulation settings
		 */

		final int bounds = 128;
		final int shape_bounds = 64;
		BSim sim = new BSim();		
		sim.setDt(1);//default is 0.01
		sim.setTimeFormat("0.00");
		sim.setSolid(true,true,true);
		sim.setBound(bounds,bounds,bounds);
		sim.setSimulationTime(400);		// Total simulation time [sec]




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
		final double firstwave = 900; // Time slot where first wave stops
		final double secondwave = 1440; // Time slot for second wave beginning
		final BSimChemicalField field = new BSimChemicalField(sim, new int[]
				{xyboxno,xyboxno,boxno}, diffusivity, decayRate);
//		field.linearZ(0,c);
		final double width = 0.4;

		
		//Set temporary number and total number of bacteria
		
		long totbac = Math.round(11459 - 2.1972*(Math.pow(width*shape_bounds, 2)));;
		
		//For even space
		double tempnum = 1/(2*totbac);
		//For packed
		//double tempnum = 1;
		
		/*********************************************************
		 * Set up the bacteria
		 */
		
		
		
		final Vector<BSimBacterium> bacteria = new Vector<BSimBacterium>();		
		while(bacteria.size() < totbac) {
			
		
			////////////////////////////////////////////////////////////////////////////////
			//////////////////// Generate position of bacteria /////////////////////////////
			////////////////////////////////////////////////////////////////////////////////
			 double hole_size = 0;
			 double u = Math.random();
			 double theta = 2 * Math.PI * u;

			
			 // Full box
			 double x0 = Math.random();
			 double y0 = Math.random();
			 double z0 = Math.random();
			 
			 while ( 0.5-(width/2) < x0 && x0 < 0.5+(width/2) && 0.5-(width/2) < z0  && z0 < 0.5+(width/2)) {
				 x0 = Math.random();
				 y0 = Math.random();
				 z0 = Math.random();
			 }

			 // Scale the shape to fit the box
			 
			 // double rand = 0.2;
			 double x00 = (x0 * (1) * shape_bounds); 		 
			 double y00 =  y0 * (1) * shape_bounds;
			 double z00 =  z0 * (1) * shape_bounds;
			 double rad = 0.05;
			 
			 
			////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////
			 
			
			BSimBacterium p = new BSimBacterium(sim, 
				//This uses random bacteria	

					//Spherical Test

					new Vector3d(0.25 * bounds+(x00), 
							0.25 * bounds + y00, 
							0.25 * bounds + z00)) {
				
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
					if ((sim.getTime()%1440 > 0 &&
						 sim.getTime()%1440 < 4 &&
						 sim.getTime() < secondwave))
						{
						field.addQuantity(bounds/4,3*bounds/4,bounds/4, quantityadd);
											
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
		
		String resultsDir = BSimUtils.generateDirectoryPath("./results/2D/newdata/Square hole/");
		
		 // Name the file
			BSimLogger logger = new BSimLogger(sim, resultsDir + "SquareHoleTest.csv") {
			
			@Override
			
			public void before() {
				super.before();
				// Write a header containing the names of variables we will be exporting
				//write("%time,%Intensity"); 
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
                
                
                // Define an array of strings to hold an array of strings
                String[] intensitiesPlane = new String[shape_bounds];
                
                // Initialise all the strings in the array as empty strings
                for (int index = 0; index < shape_bounds; index++)
                	intensitiesPlane[index] = "";
                


                
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



                }
                else{ // give the intensities along the z-axis
                		// Iterate over all the boxes in the x dimension
                    	for (int Xval = 0; Xval < 64; Xval++ ) {
                    		
                			for (int Zval = 0; Zval < 64; Zval++) {
                				//Intensity = field.getConc(0, 0, i);
                				totchem += "," + field.totalQuantity();
                				// adjust the position of the boxes to frame the circle
                				// The sphere is centred in the middle of the box
                				intensities += "," + (field.getConc(32 + Xval, 96, 32+Zval)); // these are supposed to be the indexes of the boxes coordinates, ie boxsize dependent
                			}
                			
                			// add the string of y intensities to the array of x values
                			intensitiesPlane[Xval] = intensities;
                			// remove all previous values
                			intensities = "";
                    	}

                }
                
                /*
                 * Put all the data in one line to be processed by Python later
                 * Allows each time step to be on its own line.
                 * */
                
                
                for (int index1 = 0; index1 < intensitiesPlane.length; index1++) {
                	//totalIntensities = totalIntensities + intensitiesPlane[index1];
                	write(sim.getFormattedTime() + intensitiesPlane[index1]);
                }


                //}

                // Write the time and number of collisions to file

            }
        };
        sim.addExporter(logger);
        

        // Run the simulation

        //sim.export();
        sim.preview();


    }
}
