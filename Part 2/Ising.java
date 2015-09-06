import java.awt.*;
import java.util.Scanner;
import java.awt.event.*;
import java.text.DecimalFormat;
import java.lang.Math.*;
import java.io.*;

class Ising extends Canvas implements Runnable {
    
    int size = 100;                             // Initial No. point per lattice side
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
    Button startButton = new Button("  Start  ");
    Button stopButton = new Button ("  Stop   ");
    Button dynamicsButton = new Button("Dynamics: Single");
    Scrollbar p1Scroller;                        // scrollbar to adjust temperature
    Scrollbar sizeScroller;                     // scrollbar to adjust size
    Scrollbar p2Scroller;
    Scrollbar p3Scroller;
    
    Label tLabel = new Label("P1 = 0.5" + "       ");	// text label next to scrollbar
    Label sizeLabel = new Label ("Size = " + size + "x" + size ) ;
    DecimalFormat twoPlaces = new DecimalFormat("0.00");	// to format temperature readout
    Image offScreenImage;                       // for double-buffering
    Graphics offScreenGraphics;
    Color outOfBounds = new Color(0,0,0,0); // Transp.
    
    // Constructor method handles all the initializations:
    // Set up the GUI...
    Ising() {
        
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
        
        controlPanel.add(tLabel);
        controlPanel.add(new Label("     "));
        
        //p1 Scrollbar & handlers (0.01K < T < 4.0K)
        p1Scroller = new Scrollbar(Scrollbar.HORIZONTAL,50,1,0,101) {
            public Dimension getPreferredSize() {
                return new Dimension(100,15);
            }
        };
        p1Scroller.setBlockIncrement(1);
        p1Scroller.addAdjustmentListener(new AdjustmentListener() {
            public void adjustmentValueChanged(AdjustmentEvent e) {
                tLabel.setText("P1 = " + twoPlaces.format(p1Scroller.getValue()/100.0));
            }
        });
        
        controlPanel.add(p1Scroller);
        controlPanel.add(new Label("     "));
        
        //p2 Scrollbar & handlers ( 0 < p2 < 1 )
        p2Scroller = new Scrollbar(Scrollbar.HORIZONTAL,50,1,0,101) {
            public Dimension getPreferredSize() {
                return new Dimension(100,15);
            }
        };
        p2Scroller.setBlockIncrement(1);
        p2Scroller.addAdjustmentListener(new AdjustmentListener() {
            public void adjustmentValueChanged(AdjustmentEvent e) {
                tLabel.setText("P2 = " + twoPlaces.format(p2Scroller.getValue()/100.0));
            }
        });
        
        controlPanel.add(p2Scroller);
        controlPanel.add(new Label("     "));
        
        //p2 Scrollbar & handlers ( 0 < p2 < 1 )
        p3Scroller = new Scrollbar(Scrollbar.HORIZONTAL,50,1,0,101) {
            public Dimension getPreferredSize() {
                return new Dimension(100,15);
            }
        };
        p3Scroller.setBlockIncrement(1);
        p3Scroller.addAdjustmentListener(new AdjustmentListener() {
            public void adjustmentValueChanged(AdjustmentEvent e) {
                tLabel.setText("P3 = " + twoPlaces.format(p3Scroller.getValue()/100.0));
            }
        });
        
        controlPanel.add(p3Scroller);
        controlPanel.add(new Label("     "));
        
        //SizeScroller & handlers (20 < Rows/Columns < 200)
        sizeScroller = new Scrollbar(Scrollbar.HORIZONTAL,101,1,20,101) {
            public Dimension getPreferredSize() {
                return new Dimension(100,15);
            }
        };
        
        controlPanel.add(sizeLabel);
        sizeLabel.setAlignment(Label.RIGHT);
        
        sizeScroller.setBlockIncrement(1);
        
        sizeScroller.addAdjustmentListener(new AdjustmentListener() {
            public void adjustmentValueChanged(AdjustmentEvent e) {
                sizeLabel.setText("Size = " + sizeScroller.getValue() + " X " + sizeScroller.getValue());
                size = sizeScroller.getValue();
                N = size*size;
                
                offScreenImage = createImage(canvasSize,canvasSize);
                offScreenGraphics = offScreenImage.getGraphics();
                
                // initialize the lattice
                for (int i=0; i<size; i++) {
                    for (int j=0; j<size; j++) {
                        if (Math.random() < 0.9) s[i][j] = 1; else s[i][j] = -1;
                        colorSquare(i,j);
                    }
                }
                
                // For loop when changing size of lattice
                for (int i=size+1; i<101; i++) {
                    for (int j=size+1; j<101; j++) {
                        killColor(i,j);
                    }
                }
                
            }
        });
        
        controlPanel.add(sizeScroller);
        controlPanel.add(new Label("     "));
        
        
        //Add buttons to interface (start/pause/resume, dynamics type & stop)
        startButton.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                running = !running;
                if (running) startButton.setLabel("Pause"); else startButton.setLabel("Resume");
                stopButton.setLabel("Stop");
                sizeScroller.setVisible(false);
                
            }
        });
        controlPanel.add(startButton);
        
        stopButton.addActionListener(new ActionListener(){
            public void actionPerformed(ActionEvent e) {
                
                
                if (running){
                    running = !running;
                    stopButton.setLabel("Change size and run");
                    sizeScroller.setVisible(true);
                }
                else {
                    
                    offScreenImage = createImage(canvasSize,canvasSize);
                    offScreenGraphics = offScreenImage.getGraphics();
                    
                    sizeScroller.setVisible(false);
                    
                    for (int i=0; i<size; i++) {
                        for (int j=0; j<size; j++) {
                            if (Math.random() < 0.99) s[i][j] = 1; else s[i][j] = -1;
                            colorSquare(i,j);
                        }
                    }
                    startButton.setLabel("Start");
                }
            }
        });
        controlPanel.add(stopButton);
        
        dynamicsButton.addActionListener(new ActionListener(){
            public void actionPerformed(ActionEvent e) {
                
                if (!kawasakiDynamics)
                    dynamicsButton.setLabel("Dyn: Lattice");
                else{
                    dynamicsButton.setLabel("Dyn: Single");
                }
                // Change dynamics type
                kawasakiDynamics = !kawasakiDynamics;
            }
        });
        controlPanel.add(dynamicsButton);
        
        isingFrame.pack();
        offScreenImage = createImage(canvasSize,canvasSize);
        offScreenGraphics = offScreenImage.getGraphics();
        
        for (int i=0; i<size; i++) {					// initialize the lattice the first time
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
            
            if (running) {
                double p1 = p1Scroller.getValue() / 100.0;  // Chance of getting infected              (S -> I)
                double p2 = p2Scroller.getValue() / 100.0; // Chance of getting recovered             (I -> R)
                double p3 = p3Scroller.getValue() / 100.0; // Chance of getting susceptible again     (R -> S)
                
                // if Lattice Updating
                if (kawasakiDynamics){
                    
                    for (int i=0; i<size; i++) {					// initialize the lattice the first time
                        for (int j=0; j<size; j++) {
                            
                            int count = 0;
                            
                            int above, below;
                            int left, right;
                            
                            if (i == 0) left = size - 1; else left = i - 1;
                            if (i == size - 1 ) right = 0; else right = i + 1;
                            
                            if (j == 0) above = size - 1; else above = j - 1;
                            if (j == size - 1 ) below = 0; else below = j + 1;
                            
                            if (s[left][j]  == -1) count++;
                            if (s[right][j] == -1) count++;
                            if (s[i][above] == -1) count++;
                            if (s[i][below] == -1) count++;
                            
                            //Handle S cell
                            if (s[i][j] == 1){
                                
                                // S becomes I
                                if ((count > 0) && (Math.random() < p1) ){
                                    sUpdate[i][j] = -1;
                                }
                                else {
                                    sUpdate[i][j] = s[i][j];
                                }
                            }
                            
                            //Handle I cell
                            else if (s[i][j] == -1){
                                
                                if (Math.random() < p2){
                                    // I becomes R
                                    sUpdate[i][j] = 0;
                                }
                                else {
                                    sUpdate[i][j] = s[i][j];
                                }
                            }
                            
                            //Handle R cell
                            else{
                                
                                if(Math.random() < p3) {
                                    //R becomes S
                                    sUpdate[i][j] = 1;
                                }
                                else {
                                    sUpdate[i][j] = s[i][j];
                                }
                            }
                        }
                    }
                    s = sUpdate;
                    
                    updateAllColors();
                    
                    // causes update method to be called
                    repaint();
                    
                }
                // End Kawasaki
                
                // if Single Cell updating
                else {
                    
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
                    repaint();
                }
                // End Glauber
            }
            
            if (kawasakiDynamics){
                try { Thread.sleep(150); } catch (InterruptedException e) {}	// sleep time in milliseconds
            }
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
        
        new Ising();
    }
}



