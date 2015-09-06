import java.awt.*;
import java.util.Scanner;
import java.awt.event.*;
import java.text.DecimalFormat;
import java.lang.Math.*;
import java.io.*;

class Sirs extends Canvas implements Runnable {
    
    int size = 50;                             // Initial No. point per lattice side
    int squareWidth = 5;                        // pixels across one lattice site
    int canvasSize = 600;                       // total pixels across canvas
    int[][] s = new int[size][size];            // the 2D array of dipoles
    int[][] sUpdate = new int[size][size];      // the 2D array of dipoles for full lattice updating
    
    // s[][] =  1 : Susceptible :   White
    // s[][] = -1 : Infected    :   Red
    // s[][] =  0 : Recovered   :   Green
    
    int sweeps = 0;                             // #Sweeps
    int N = size*size;                          // #Cells in lattice
    int nInfected;				// # infected cells
    float psi;                                  // frac. I-cells
    float psiAve;                               // <psi>
    float psi2Ave;                              // <psi*psi>
    int count = 0;
    int cellFlip = 0;
    
    Color upColor = new Color(255,255,255);		// White
    Color downColor = new Color(153,223,174);	// green
    Color rColor = new Color(230, 40, 40);      // Red
    
    boolean running = false; // true when simulation is running
    boolean runAnalytics = false;
    boolean waveAnalyticsRun = false;
    boolean kawasakiDynamics = false; // false = Glauber, true = Kawasaki
    
    Label p1Label = new Label(" P1: 0.0 ");
    Label p2Label = new Label (" P2: 0.5 ");
    Label p3Label = new Label(" P3: 0.0 ");
	Label progress = new Label( "" );

	DecimalFormat twoPlaces = new DecimalFormat("0.0");

    Image offScreenImage;                       // for double-buffering
    Graphics offScreenGraphics;
    Color outOfBounds = new Color(0,0,0,0); // Transp.
    
    // Constructor method handles all the initializations:
    // Set up the GUI...
    Sirs() {
        
        Frame isingFrame = new Frame("SIRS Model");
        isingFrame.addWindowListener(new WindowAdapter() {
            public void windowClosing(WindowEvent e) {
		    System.exit(0);            			       
            }
        });
        
        //Initialise the canvas with components
        Panel canvasPanel = new Panel();
        isingFrame.add(canvasPanel);
        canvasPanel.add(this);
        
        setSize(canvasSize,canvasSize);
        Panel controlPanel = new Panel();
        isingFrame.add(controlPanel,BorderLayout.SOUTH);
        
        controlPanel.add(p1Label);
        controlPanel.add(p2Label);
	controlPanel.add(p3Label);
        
        isingFrame.pack();
        offScreenImage = createImage(canvasSize,canvasSize);
        offScreenGraphics = offScreenImage.getGraphics();
        
	// initialize the lattice the first time
        for (int i=0; i<size; i++) {					
            for (int j=0; j<size; j++) {
                if (Math.random() < 0.99) s[i][j] = 1; else s[i][j] = -1;
                colorSquare(i,j);
            }
        }
        
        isingFrame.setVisible(true);
        
        Thread t = new Thread(this);		// create a thread to run the simulation
        t.start();
    }
    
    // Run method gets called by new thread to carry out the simulation:
    public void run() {
        while (true) {

	try{
               		
 	//Method for analysing the lattice under different p1 and p3:
	System.out.print("running analytics");

        double[] p3Data = new double[11];
	double[] p1Data = new double[11];
        double[][] contourData = new double[11][11];

	//Data for wave analysis part
	double[][] waveContourData = new double[11][11];
        
        double p1 = 0.0;
        double p2 = 0.5;
        double p3 = 0.0;
        
	
        //Sample <I> for various p1 and p3
        for (int ii=0; ii<11; ii++) { // p1
            
            p3 = 0.0;
            
            for (int jj=0; jj<11; jj++) { // p3
		p1Label.setText(" P1 = " + twoPlaces.format(p1));
		p3Label.setText(" P3 = " + twoPlaces.format(p3));

                boolean analysing = true;
                
                while(analysing){
                    
                    nInfected = 0;  
                    int count = 0;
                                     
                    int i = (int) (Math.random() * size);	// choose a random row and column
                    int j = (int) (Math.random() * size);
                    
                    int above, below;
                    int left, right;
                    
                    if (i == 0) left = size - 1; else left = i - 1;
                    if (i == size - 1 ) right = 0; else right = i + 1;
                    
                    if (j == 0) above = size - 1; else above = j - 1;
                    if (j == size - 1 ) below = 0; else below = j + 1;
                    
                    if (s[left][j] == -1) count++;
                    if (s[right][j] == -1) count++;
                    if (s[i][above] == -1) count++;
                    if (s[i][below] == -1) count++;
                    
                    //Handle S cell
                    if (s[i][j] == 1){
                        
                        // S becomes I
                        if ((count > 0) && (Math.random() < p1) ){
                            s[i][j] = -1;
                        }
                    }  
                    //Handle I cell
                    else if (s[i][j] == -1) {
                        if (Math.random() < p2){
                        // I becomes R
                        s[i][j] = 0;
			}
                    }
                    
                    //Handle R cell
                    else{
                        if(Math.random() < p3){
                        //R becomes S
                        s[i][j] = 1;
			}
                    }
                    
                    colorSquare(i,j);
                    // causes update method to be called
                                        
                    for (int xx=0; xx<size; xx++) {
                        for (int yy=0; yy<size; yy++) {
                            if (s[xx][yy] == -1) {
                                nInfected ++;
                            }
                        }
                    }
                    
                    cellFlip ++;
                    psi += (float)nInfected / (float)N;

		// One MC-Sweep: cellFlip == N
		if (cellFlip == N){
                       
			//At this point psi = # i-cells for the last MC-sweep
                        sweeps++;

			// Using average psiAve of last 100 as the 'actual' value: To ensure system is in stability region
			// Exit current analysing cycle
                        if ((sweeps > 200) || (nInfected == 0)){

                            psiAve = psiAve + psi; // summing i-cells
			    psi2Ave = psi2Ave + psi*psi; // summing (i-cells)**2

			if ((sweeps == 300) || (nInfected == 0)){
				if (sweeps == 300){
                               		psiAve = psiAve / (float)(100);
					psi2Ave = psi2Ave / (float)(100);
                       		 }
				if (nInfected == 0){
					psiAve = 0;
					psi2Ave = 0;
				}
			// Handle cycle-exit: For absorbing-state systems or systems in dynamic equil.: Save data to arrays
							
		            analysing = false;
                            contourData[ii][jj] = psiAve / (float)N;
			    waveContourData[ii][jj] = (psi2Ave - psiAve*psiAve) / (float)N ;
                            System.out.println("(" + p1 + ", " + p3 + ")" + " :" + contourData[ii][jj] + " : " + waveContourData[ii][jj]);

			 // initialize the new lattice
                            for (int ig=0; ig<size; ig++) {
                                for (int jg=0; jg<size; jg++) {
                                    if (Math.random() < 0.99) s[ig][jg] = 1; else s[ig][jg] = -1;
                                    colorSquare(ig,jg);
                                }
                           }

			//Reset Sweeps etc.
			sweeps = 0;
                        psiAve = 0;
			psi2Ave = 0;
	}			
} // End Exit cycle

			// Else new MC Cycle
			psi = 0;	
                        cellFlip = 0;

                        }// End MC-iterator

                    repaint();

                }//end analysing while-loop
                p3 += 0.1;
            }
            p1 += 0.1;
        }   

	//When for-loops have finished, print data to file:
	PrintWriter output = new PrintWriter(new FileWriter("data.txt"));
	PrintWriter waveOutput = new PrintWriter(new FileWriter("waveData.txt"));

	for (int lineX = 0; lineX < 11; lineX++){
	    for(int lineY = 0; lineY < 11; lineY++){

		output.print(contourData[lineX][lineY] + " ");
		waveOutput.print(waveContourData[lineX][lineY] + " ");


	    }
	    output.println();
	    waveOutput.println();
	}

	// Closes the outputfile. //
	output.close();
	waveOutput.close();	

		}
			catch(IOException ex){    
		}

System.exit(0);
    }         
}
           
    
    
    //For full update scheme
    void updateAllColors() {
        
        for (int i=0; i<size; i++) {
            for (int j=0; j<size; j++) {
                colorSquare(i,j);
            }
        }
    }
    
    // Color a given square according to the site's orientation:
    void colorSquare(int i, int j) {
        if (s[i][j] == 1) offScreenGraphics.setColor(upColor);
        else if (s[i][j] == 0) offScreenGraphics.setColor(downColor);
        else offScreenGraphics.setColor(rColor);
        offScreenGraphics.fillRect(i*squareWidth,j*squareWidth,squareWidth,squareWidth);
    }
    
    //Transparetn Color (for out of bounds cells)
    void killColor(int i, int j){
        offScreenGraphics.setColor(outOfBounds);
        offScreenGraphics.fillRect(i*squareWidth,j*squareWidth,squareWidth,squareWidth);
        
    }
    
    // Override default update method to skip drawing the background:
    public void update(Graphics g) {
        paint(g);
    }
    
    // Paint method just blasts the off-screen image to the screen:
    public void paint(Graphics g) {
        g.drawImage(offScreenImage,0,0,canvasSize,canvasSize,this);
    }
    
    // Main method just calls constructor to get started:
    public static void main(String[] args) throws IOException {
        
        new Sirs();
    }
}



