import java.awt.*;
import java.util.Scanner;
import java.awt.event.*;
import java.text.DecimalFormat;
import java.lang.Math.*;
import java.io.*;

class Waves extends Canvas implements Runnable {
    
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
    Label p1Label = new Label(" Sweep: ");
    
    Image offScreenImage;                       // for double-buffering
    Graphics offScreenGraphics;
    Color outOfBounds = new Color(0,0,0,0); // Transp.
    
    // Constructor method handles all the initializations:
    // Set up the GUI...
    Waves() {
        
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
                
                float[] waveStream = new float[501];
                int[] sweepCounter  = new int[501];
                
                //Wave setup
                double p1 = 0.8;
                double p2 = 0.1;
                double p3 = 0.01;
                
                while( sweeps < 500 ){
                    
                    p1Label.setText(" Sweep: " + (sweeps+1) + " / 5000.");
                    
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
                        
                        waveStream[sweeps] = (float)psi / (float)N;
                        sweepCounter[sweeps] = sweeps ;
                        
                        // Else new MC Cycle
                        psi = 0;
                        cellFlip = 0;
                        
                    }// End MC-iterator
                    
                    repaint();
                    
                } // End for-loop: Write data to file
                
                PrintWriter wavefunc = new PrintWriter(new FileWriter("waveFunc.txt"));
                for(int line = 0; line < 501; line++){
                    wavefunc.print(sweepCounter[line] + " " +  waveStream[line]);
                    wavefunc.println();
                }
                wavefunc.close();
                
                
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
        
        new Waves();
    }
}



