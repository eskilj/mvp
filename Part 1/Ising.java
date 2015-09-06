import java.awt.*;
import java.util.Scanner;
import java.awt.event.*;
import java.text.DecimalFormat;
import java.lang.Math;

class Ising extends Canvas implements Runnable {
    
//    
//    // Set size of the 2D lattice
//    Scanner inputInt = new Scanner(System.in);
//    int size = inputInt.nextInt();              // number of lattice sites in a row (change if desired)
    
    int size = 100;                             // Initial No. point per lattice side
	int squareWidth = 5;                        // pixels across one lattice site
    int arraySize = 100;
	int canvasSize = 600;                       // total pixels across canvas
    int[][] s = new int[size][size];            // the 2D array of dipoles (each equal to 1 or -1)

    int nData = 0;                              // Data counter
    int n = 0;                                  // Total Monte-Carlo Counter
    
    int M;
    double absSumM;                             // <|M|>
    double sumMa;                               // <M>
    double sumMaMa;                             // <M*M>
    
    int E;
    double sumEn;                               // <E>
    double sumEnEn;                             // <E*E>
    
    double c;                                   // Heat Capacity
    double chi;                                 // Chi
    
    int[] sumM = new int[arraySize];            // Array used for graphing
    int[] sumE = new int[arraySize];            // Array used for graphing

    boolean running = false; // true when simulation is running
    boolean kawasakiDynamics = false; // false = Glauber, true = Kawasaki
	Button startButton = new Button("  Start  ");
    Button stopButton = new Button ("  Stop   ");
    Button dynamicsButton = new Button("Dyn: Glauber");
	Scrollbar tScroller;                        // scrollbar to adjust temperature
    Scrollbar sizeScroller;                     // scrollbar to adjust size
    
	Label tLabel = new Label("Temperature = 2.27");	// text label next to scrollbar
    Label sizeLabel = new Label ("Size = " + size + "x" + size ) ;
	DecimalFormat twoPlaces = new DecimalFormat("0.00");	// to format temperature readout
	Image offScreenImage;                       // for double-buffering
	Graphics offScreenGraphics;
	Color upColor = new Color(0,0,0);		// black
	Color downColor = new Color(255,255,255);	// white
    Color outOfBounds = new Color(0,0,0,0); // Transp.
	
	// Constructor method handles all the initializations:
    // Set up the GUI...
	Ising() {
        
		Frame isingFrame = new Frame("Ising Model");
        isingFrame.add(this);
        setSize(900, 600);
		isingFrame.addWindowListener(new WindowAdapter() {
			public void windowClosing(WindowEvent e) {
				System.exit(0);							// close button exits program
			}
		});
        
        //Initialise the canvas with components
		Panel canvasPanel = new Panel();
		isingFrame.add(canvasPanel, BorderLayout.CENTER);
		canvasPanel.add(this);
		setSize(canvasSize,canvasSize);
        
		Panel controlPanel = new Panel();
        controlPanel.setLayout(new GridLayout(20,0));
		isingFrame.add(controlPanel,BorderLayout.EAST);

        
		controlPanel.add(tLabel);
        controlPanel.add(new Label(""));
        
        //Temperature Scrollbar & handlers (0.01K < T < 4.0K)
		tScroller = new Scrollbar(Scrollbar.HORIZONTAL,227,1,1,401) {
			public Dimension getPreferredSize() {
				return new Dimension(75,15);
			}
		};
		tScroller.setBlockIncrement(1);
		tScroller.addAdjustmentListener(new AdjustmentListener() {
			public void adjustmentValueChanged(AdjustmentEvent e) {
				tLabel.setText("Temp = " + twoPlaces.format(tScroller.getValue()/100.0));
			}
		});
        
		controlPanel.add(tScroller);
		controlPanel.add(new Label(""));
        
        //SizeScroller & handlers (20 < Rows/Columns < 200)
        sizeScroller = new Scrollbar(Scrollbar.HORIZONTAL,200,1,20,201) {
            public Dimension getPreferredSize() {
                return new Dimension(75,15);
            }
        };
        
        controlPanel.add(sizeLabel);
        sizeLabel.setAlignment(Label.RIGHT);
        
        sizeScroller.setBlockIncrement(1);

        sizeScroller.addAdjustmentListener(new AdjustmentListener() {
        public void adjustmentValueChanged(AdjustmentEvent e) {
        sizeLabel.setText("Size = " + sizeScroller.getValue() + " X " + sizeScroller.getValue());
            
            size = sizeScroller.getValue();
            
            offScreenImage = createImage(canvasSize,canvasSize);
            offScreenGraphics = offScreenImage.getGraphics();
            
            // initialize the lattice
            for (int i=0; i<size; i++) {
                for (int j=0; j<size; j++) {
                    if (Math.random() < 0.5) s[i][j] = 1; else s[i][j] = -1;
                    colorSquare(i,j);
                }
            }
            
            // For loop when changing size of lattice
            for (int i=size+1; i<201; i++) {
                for (int j=size+1; j<201; j++) {
                    killColor(i,j);
                }
            }
         
        }
        });
        
        controlPanel.add(sizeScroller);
        controlPanel.add(new Label(""));
        
        
        //Add buttons to interface (start/pause/resume, dynamics type & stop)
		startButton.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				running = !running;
				if (running) startButton.setLabel("Pause"); else startButton.setLabel("Resume");
                sizeScroller.setVisible(false);
			}
		});
		controlPanel.add(startButton);
        
        stopButton.addActionListener(new ActionListener(){
            public void actionPerformed(ActionEvent e) {
            
                if (running) running = !running;
                startButton.setLabel("Start");
                sizeScroller.setVisible(true);
            
            }
        });
        controlPanel.add(stopButton);
        
        dynamicsButton.addActionListener(new ActionListener(){
            public void actionPerformed(ActionEvent e) {
                
                if (!kawasakiDynamics)
                dynamicsButton.setLabel("Dyn: Kawasaki");
                else{
                dynamicsButton.setLabel("Dyn: Glauber");
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
				if (Math.random() < 0.5) s[i][j] = 1; else s[i][j] = -1;
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
				double temp = tScroller.getValue() / 100.0;
				for (int step=0; step<1000; step++) {		// number of steps (1000)
                    
                    int N = size*size;
                    
                    // if kawasaki
                    if (kawasakiDynamics){
                    
                        // First site
                        int i1 = (int) (Math.random() * size);
                        int j1 = (int) (Math.random() * size);
                        int siteSpin1 = s[i1][j1]; // current site config
                        
                        //Second site
                        int i2 = (int) (Math.random() * size);
                        int j2 = (int) (Math.random() * size);
                        int siteSpin2 = s[i2][j2];
                        
                        if ((s[i1][j1] != s[i2][j2])){ // current site config
                            double eDiffKawasaki = deltaU(i1,j1) + deltaU(i2, j2);
                            if ((eDiffKawasaki <= 0) || (Math.random() < Math.exp(-eDiffKawasaki/temp))) {
                                
                                s[i1][j1] = siteSpin2;
                                colorSquare(i1,j1);
                                
                                s[i2][j2] = siteSpin1;
                                colorSquare(i2,j2);
                            }
                        }
                    }
                    // End Kawasaki
                    
                    // if Glauber
                    else {
                    
                        int i = (int) (Math.random() * size);	// choose a random row and column
                        int j = (int) (Math.random() * size);
                        double eDiff = deltaU(i,j);				// compute energy change if flipped
                        if ((eDiff <= 0) || (Math.random() < Math.exp(-eDiff/temp))) {
                            s[i][j] *= -1;
                            colorSquare(i,j);
                        }
                        

                    }
                    // End Glauber
                    
                    // Total energy E
                    E = 0;
                    
                    for (int x = 0; x < size; x++) {
                        int right = x + 1;
                        if (right == size)
                            right = 0;
                        for (int y = 0; y < size; y++) {
                            int up = y + 1;
                            if (up == size)
                                up = 0;
                            int sum = s[x][up] + s[right][y];
                            E += - s[x][y] * sum;
                        }
                    }

                    // Total Magnitisation
                    M = 0;
                    for(int x = 0; x < size; x++){
                        for(int y = 0;y < size; y++){
                            
                            M += s[x][y];
                        }
                    }
                    
                    double norm = 1.0 / N;
                    
                    if (nData > 0 ){
                    
                        norm /= nData;
                    }
                    
                    absSumM += Math.abs(M);
                    sumMa += M;
                    sumMaMa += M * (double)M;
                    
                    sumEn += E;
                    sumEnEn += E * E;
                    nData++;
                    
                    c = ((sumEnEn*norm) - (sumEn*norm*sumEn*norm)) / (temp * temp * N);
                    chi = (sumMaMa*norm - N*sumMa*sumMa*norm*norm) / (temp);
                    
//                    System.out.println("cv = " + c);
//                    System.out.println("chi = " + chi);
                    
                    n++;
                    
				}
                // causes update method to be called
				repaint();
                
 
			}
			try { Thread.sleep(1); } catch (InterruptedException e) {}	// sleep time in milliseconds
		}
	}

    // Given a lattice site, compute hypothetical energy change from flip
    double deltaU(int i, int j) {
        int leftS = 0, rightS = 0, topS = 0, bottomS = 0;  // values of neighboring spins
        if (i == 0) leftS = s[size-1][j]; else leftS = s[i-1][j];
        if (i == size-1) rightS = s[0][j]; else rightS = s[i+1][j];
        if (j == 0) topS = s[i][size-1]; else topS = s[i][j-1];
        if (j == size-1) bottomS = s[i][0]; else bottomS = s[i][j+1];
        return 2.0 * s[i][j] * (leftS + rightS + topS + bottomS);
    }

	// Color a given square according to the site's orientation:
	void colorSquare(int i, int j) {
		if (s[i][j] == 1) offScreenGraphics.setColor(upColor); 
					 else offScreenGraphics.setColor(downColor);
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
    public static void main(String[] args) {
        
		new Ising();    
    }  
}



