import java.applet.*;
import java.awt.*;
import java.awt.event.*;
import java.awt.geom.*;
import java.awt.image.*;
import static java.lang.Math.*;
import java.text.*;
import javax.swing.*;
import java.util.Hashtable;
import java.util.Random;
import java.io.*;

public class LJ3MDApp extends JFrame implements ActionListener
{
    double rhoInit = 0.5;
    LJ3MD md = new LJ3MD(rhoInit);
    int delay = 100;
    int speed = 10; // md steps per frame
    Timer timer;
    
    XYZCanvas canvas; // 3D drawing here
    JPanel cpnl, spnl;
    JTextField tNum       = new JTextField("" + md.N);
    JTextField tTemp      = new JTextField("1.0");
    JTextField tRho       = new JTextField("" + rhoInit);
    JTextField tSpeed     = new JTextField("" + speed);
    JButton    bReset     = new JButton("Reset");
    JButton    bRetime    = new JButton("Reset time");
    JToggleButton bStart  = new JToggleButton("Start");
    JCheckBox  bTstat     = new JCheckBox("Thermostat", true);
    JCheckBox  bPot       = new JCheckBox("L-J Potential", true);
    JLabel     lStatus    = new JLabel("Status");
    JTextField tAvK       = new JTextField("0");
    JTextField tAvU       = new JTextField("0");
    JTextField tAvp       = new JTextField("0");
    public static final long serialVersionUID = 1L;
    
    public LJ3MDApp() {
        tNum.setHorizontalAlignment(JTextField.CENTER);
        tTemp.setHorizontalAlignment(JTextField.CENTER);
        tRho.setHorizontalAlignment(JTextField.CENTER);
        tSpeed.setHorizontalAlignment(JTextField.CENTER);
        
        tAvK.setHorizontalAlignment(JTextField.RIGHT);
        tAvU.setHorizontalAlignment(JTextField.RIGHT);
        tAvp.setHorizontalAlignment(JTextField.RIGHT);
        
        float[] aveKing = new float[501];
        float[] avePot = new float[501];
        float[] aveEn = new float[501];
        
        
        JFrame box = new JFrame();
        box.setLayout(new BorderLayout());
        box.setSize(1000, 1000);
        box.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        
        cpnl = new JPanel(); // create a panel for controls
        cpnl.setLayout(new GridLayout(18, 2));
        box.add(cpnl, BorderLayout.EAST);
        
        // add controls
        cpnl.add(bStart);
        bStart.addActionListener(this);
        
        cpnl.add(bReset);
        bReset.addActionListener(this);
        
        cpnl.add(new JLabel(" N:"));
        tNum.addActionListener(this);
        cpnl.add(tNum);
        
        cpnl.add(new JLabel(" Density (\u03c1):"));
        tRho.addActionListener(this);
        cpnl.add(tRho);
        
        cpnl.add(new JLabel(" Steps/frame:"));
        tSpeed.addActionListener(this);
        cpnl.add(tSpeed);
        
        cpnl.add(bTstat);
        bTstat.addActionListener(this);
        
        cpnl.add(bPot);
        bPot.addActionListener(this);
        
        cpnl.add(new JLabel(" < K/N > :"));
        tAvK.setEditable(false);
        cpnl.add(tAvK);
        
        cpnl.add(new JLabel(" Temperature:"));
        tTemp.setEditable(false);
        cpnl.add(tTemp);
        
        cpnl.add(new JLabel(" < U/N > :"));
        tAvU.setEditable(false);
        cpnl.add(tAvU);
        
        cpnl.add(new JLabel(" < pressure > :"));
        tAvp.setEditable(false);
        cpnl.add(tAvp);
        
        cpnl.add(bRetime);
        bRetime.addActionListener(this);
        
        spnl = new JPanel(); // create a panel for status
        box.add(spnl, BorderLayout.SOUTH);
        lStatus.setFont(new Font("Courier", 0, 12));
        spnl.add(lStatus);
        
        canvas = new XYZCanvas();
        box.add(canvas, BorderLayout.CENTER);
        
        timer = new Timer(delay, this);
        timer.start();
        //        timer.stop();
        box.setVisible(true);
        
    }
    
    public void start() { }
    
    public void update(Graphics g) { paint(g); }
    
    DecimalFormat df = new DecimalFormat("###0.000");
    //    int frame = 0;
    public void paint(Graphics g) {
        //System.out.println("frame: " + (frame++));
        lStatus.setText("t = " + df.format(md.dt*md.step) + ", "
                        + "N = " + md.N + ", "
                        + "E/N = " + df.format(md.E/md.N) + ", "
                        + "U/N = " + df.format(md.U/md.N) + ", "
                        + "K/N = " + df.format(md.K/md.N) + ", "
                        + "p = " + df.format(md.p) + ";");
        tAvK.setText(df.format(md.avK.getAve()/md.N) + "  ");
        tAvU.setText(df.format(md.avU.getAve()/md.N) + "  ");
        tTemp.setText(df.format((2*md.K)/(3*(md.N - 1))) + "  ");
        tAvp.setText(df.format(md.avp.getAve()) + "  ");
        canvas.refresh(md.getXWrap(), md.N, true, false);
        cpnl.repaint();
        spnl.repaint();
        
        try{
        
            PrintWriter wavefunc = new PrintWriter(new FileOutputStream(new File("energyData.txt"), true));
            wavefunc.print(md.E/md.N + " " +  md.K/md.N + " " +  md.U/md.N );
            wavefunc.println();
        wavefunc.close();
        }
        catch(IOException ex){}
        
        try{
            
            PrintWriter tempwriter = new PrintWriter(new FileOutputStream(new File("tempData.txt"), true));
            tempwriter.print(df.format((2*md.K)/(3*(md.N - 1))));
            tempwriter.println();
            tempwriter.close();
        }
        catch(IOException ex){}
        
        
    }
    
    public void actionPerformed(ActionEvent e) {
        Object src = e.getSource();
        if (src == timer) {
            for (int i = 0; i < speed; i++) // integrate a few steps
                md.vv();
            repaint();
            return;
        }
        
        boolean adjCanvasScale = false;
        
        if (src == bTstat) md.thermostat = !md.thermostat;
        
        if (src == bPot ){
            md.ljPotential = !md.ljPotential;
            md.clearData();
            
            if (timer.isRunning()) timer.stop();
            bStart.setSelected(false);
            bStart.setText("Start");
            md.init(md.rho);
            
            
        }
        
        if (src == tTemp || src == bReset) {
            double kT = Double.parseDouble(tTemp.getText().trim());
            if (kT < 1e-8) { kT = 1e-8; tTemp.setText("  " + kT); }
            md.kT = kT;
            md.clearData();
        }
        
        if (src == tRho || src == bReset) {
            double rho = Double.parseDouble(tRho.getText().trim());
            if (rho < 1e-3) { rho = 1e-3; tRho.setText("   " + rho); }
            if (rho > 1.2) { rho = 1.2; tRho.setText("   " + rho); }
            md.setDensity(rho);
            md.clearData();
            adjCanvasScale = true;
        }
        
        if (src == tSpeed || src == bReset) {
            speed = Integer.parseInt(tSpeed.getText().trim());
            if (speed < 1) { speed = 1; tSpeed.setText("   " +speed); }
        }
        
        if (src == bRetime)
            md.clearData();
        
        if (src == bStart) {
            boolean on = bStart.isSelected();
            if (on) {
                timer.restart();
                bStart.setText("Pause");
            } else {
                timer.stop();
                bStart.setText("Resume");
            }
        }
        
        if (src == tNum) {
            int n = Integer.parseInt(tNum.getText().trim());
            if (n < 2) { n = 2; tNum.setText(" " + n); }
            md.N = n;
            md.init(md.rho);
            adjCanvasScale = true;
        }
        
        if (src == bReset) {
            if (timer.isRunning()) timer.stop();
            bStart.setSelected(false);
            bStart.setText("Start");
            md.init(md.rho);
        }
        
        canvas.refresh(md.getXWrap(), md.N, true, adjCanvasScale);
        
        repaint();
    }
    
    public static void main(String[] args){
        
        LJ3MDApp letgo = new LJ3MDApp();
        letgo.setVisible(true);
        
        
    }
}


class LJ3MD {
    static final int D = 3;
    public int N = 256;
    private int DOF = N*D - D*(D+1)/2;
    double rho = 0.256;
    double dt = 0.002; // time step for integrating Newton's equation
    public double L = 10.0; // side of the box
    public double Vol = 1000.0;
    double rc = 2.5;    // cut-off radius
    double epsilon = 0.001; //interaction strength
    public double kT = 1.0; // temperature times Boltzmann constant
    double x[][]; // position, from 0 to 1
    double v[][], f[][]; //  velocity and force
    double xwrap[][]; // wrapped coordinates
    double K, U, Us, E; // kinetic, potential, and total energy
    double Vir, p; // virial and pressure
    boolean thermostat = true; // use thermostat
    boolean ljPotential = true;
    double ushift;
    double Utail; // potential energy tail correction
    double ptail; // pressure tail correction
    int step; // simulation steps
    
    double boltK = 1.38 * Math.pow(10, -23);
    Ave avU = new Ave(), avK = new Ave(), avp = new Ave();
    
    /** constructor */
    LJ3MD(double den) { init(den); }
    
    /** initialize system */
    public void init(double den) {
        int i, j, k, id;
        
        DOF = N*D - D*(D+1)/2;
        setDensity(den);
        step = 0;
        x = new double[N][D];
        v = new double[N][D];
        f = new double[N][D];
        xwrap = new double[N][D];
        for (id = 0; id < N; id++)
            for (j = 0; j < D; j++)
                v[id][j] = 0.0;
        
        int N1 = (int) (pow(2*N, 1.0/D) + .999999); // # of particles per side (cubic lattic)
        double a = 1./N1;
        for (id = 0, i = 0; i < N1 && id < N; i++) // place particles on a cubic lattice
            for (j = 0; j < N1 && id < N; j++)
                for (k = 0; k < N1 && id < N; k++) {
                    if ((i + j + k) % 2 != 0) continue;
                    x[id][0] = a * (i + .5);
                    x[id][1] = a * (j + .5);
                    x[id][2] = a * (k + .5);
                    id++;
                }
    }
    
    /** set density and recalculate tail corrections */
    void setDensity(double den) {
        rho = den;
        Vol = N/rho;
        L = pow(Vol, 1.0/D);
        
        /* compute shifts and tail corrections */
        if ((rc = 2.5) > .5*L) rc = .5*L;
        double irc = 1./rc, irc3 = irc*irc*irc, irc6 = irc3*irc3;
        ushift = 4.0*irc6*(irc6 - 1);
        Utail = 8.0*PI/9*rho*N*(irc6 - 3)*irc3;
        ptail = 16.0*PI/9*rho*rho*(irc6 - 1.5)*irc3;
        
        clearData();
    }
    
    void clearData() {
        step = 0;
        avU.clear();
        avK.clear();
        avp.clear();
    }
    
    /** integrate Newton's equation by one step, velocity Verlet */
    void vv() {
        double dth = dt*.5, dtl = dt/L;
        
        for (int id = 0; id < N; id++) // velocity-verlet, first half
            for (int d = 0; d < D; d++) {
                v[id][d] += f[id][d] * dth;
                x[id][d] += v[id][d] * dtl;
            }
        
        force(); // compute force
        
        for (int id = 0; id < N; id++) // velocity-verlet, second half
            for (int d = 0; d < D; d++)
                v[id][d] += f[id][d] * dth;
        
        rmcom(); // rm angular momentum however isn't good idea due to
        // pbc and minimal image convension
        K = vrescale(0.02); // compute kinetic energy and adjust velocity
        E = K + Us;
        step++;
        avU.add(U);
        avK.add(K);
        avp.add(p);
    }
    
    /** calculate potential energy and force
     *  half-box cutoff, minimal distance image */
    double xij[] = new double[D];
    void force() {
        int i, j, d, prcnt = 0;
        double invr2, invr6, r2, fscal, tmp;
        
        // clear the force and energy
        U = 0.0;
        Vir = 0.0;
        for (i = 0; i < N; i++)
            for (d = 0; d < D; d++)
                f[i][d] = 0.0;
        
        // loop over pairs of particles and compute energy (1) Lennard-Jones, (2) WCA;
        
        // LJ 12-6
        
            for (i = 0; i < N - 1; i++) {
                for (j = i + 1; j < N; j++) {
                    rv3Diff(xij, x[i], x[j]);
                    vpbc(xij);
                    r2 = vsqr(xij);
                    
                    //(1)
                    if (ljPotential){
                    
                        if (r2 > rc*rc) continue; // if r > rc: break
                        invr2 = 1.0/r2;
                        invr6 = invr2 * invr2 * invr2;
                        fscal = invr6  *(48.0 * invr6 - 24.0);
                        Vir += fscal;
                        fscal *= invr2;
                        for (d = 0; d < D; d++) {
                            tmp = xij[d] * fscal;
                            f[i][d] += tmp;
                            f[j][d] -= tmp;
                        }
                        U += 4.0  * invr6 * (invr6 - 1.0);
                    
                    }
                    
                    //(2)
                    else{
                        
                        if (r2 > 1.2599210492*rc*rc) {
                        
                        }
                        
                        else {
                        invr2 = 1.0/r2;
                        invr6 = invr2 * invr2 * invr2;
                        fscal = invr6  *(48.0 * invr6 - 24.0);
                        Vir += fscal;
                        fscal *= invr2;
                        for (d = 0; d < D; d++) {
                            tmp = xij[d] * fscal;
                            f[i][d] += tmp;
                            f[j][d] -= tmp;
                        }
                        U += 4.0*invr6*(invr6 - 1.0) + epsilon;
                        }
                    }
                    
                    prcnt++;
                }
            }
        
        
        Us = U - prcnt*ushift;
        U += Utail;
        p = rho*kT + Vir/(3*Vol) + ptail;
        
        
        
    }
    
    Random rng = new Random();
    
    /** velocity rescaling thermostat */
    double vrescale(double dt) {
        double ek1 = getEKin();
        
        
        
        if (!thermostat) return ek1;
        double ekav = .5f*kT*DOF;
        double amp = 2*sqrt(ek1*ekav*dt/DOF);
        double ek2 = ek1 + (ekav - ek1)*dt + amp*rng.nextGaussian();
        if (ek2 < 0) ek2 = ek1;
        double s = sqrt(ek2/ek1);
        for (int i = 0; i < N; i++)
            for (int d = 0; d < D; d++)
                v[i][d] *= s;
        return ek2;
    }
    
    /** remove motion of the center of mass */
    void rmcom() {
        for (int d = 0; d < D; d++) { // remove the center of mass motion
            double vc = 0.0;
            for (int id = 0; id < N; id++)
                vc += v[id][d];
            vc /= N;
            for (int id = 0; id < N; id++)
                v[id][d] -= vc;
        }
    }
    
    private void rv3Diff(double c[], double a[], double b[]) {
        c[0] = a[0] - b[0];
        c[1] = a[1] - b[1];
        c[2] = a[2] - b[2];
    }
    
    private double rv3Dot(double a[], double b[]) {
        return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
    }
    
    /** return the kinetic energy */
    double getEKin() {
        double ek = 0.0;
        for (int i = 0; i < N; i++)
            ek += rv3Dot(v[i], v[i]);
        return .5*ek;
    }
    
    /** image with the minimal distance */
    double pbc(double a) { return L * (a - ((int)(a + 1000.5) - 1000)); }
    void vpbc(double v[]) { v[0] = pbc(v[0]); v[1] = pbc(v[1]); v[2] = pbc(v[2]); }
    double vsqr(double v[]) { return v[0]*v[0]+v[1]*v[1]+v[2]*v[2]; }
    
    /** wrap coordinates back into box */
    public double [][] getXWrap() {
        for (int i = 0; i < N; i++)
            for (int d = 0; d < D; d++) {
                double xx = x[i][d] + 1000.0;
                xwrap[i][d] = (xx - (int) xx) * L;
            }
        return xwrap;
    }
}




/** Average and standard deviation */
class Ave {
    double cnt, xsm, x2sm;
    
    public void clear() {
        cnt = xsm = x2sm = 0;
    }
    
    public void add(double x) {
        cnt += 1;
        xsm += x;
        x2sm += x*x;
    }
    
    public double getAve() {
        return (cnt > 0) ? xsm/cnt : 0;
    }
    
    public double getVar() {
        if (cnt <= 1) return 0;
        double av = xsm/cnt;
        return x2sm/cnt - av*av;
    }
    
    public double getStd() {
        return Math.sqrt( getVar() );
    }
}

/** Panel for drawing particles */
class XYZCanvas extends JPanel
implements MouseListener, MouseMotionListener, MouseWheelListener {
    Image img;
    Graphics imgG;
    Dimension imgSize;
    
    double realSpan; // size of a real span (saved as a copy)
    double zoomScale = 1.0;
    public Matrix3D viewMatrix = new Matrix3D(); // view matrix
    private Matrix3D tmpMatrix = new Matrix3D(); // temporary matrix
    int mouseX, mouseY; // mouse position
    XYZModel model;
    
    public static final long serialVersionUID = 2L;
    
    public XYZCanvas() {
        super();
        addMouseListener(this);
        addMouseMotionListener(this);
        addMouseWheelListener(this);
    }
    
    /** Prepare a buffer for the image, return if the canvas is ready
     *  Only need to call this when the size of the canvas is changed,
     *    Since we automatically detect the size change in paint()
     *    only call this on start */
    public boolean newImgBuf() {
        Dimension sz = getSize();
        if (sz.width == 0 || sz.height == 0)
            return false;
        // quit if the current image already has the right size
        if (img != null && imgG != null && sz.equals(imgSize))
            return true;
        img = createImage(sz.width, sz.height);
        if (imgG != null) imgG.dispose();
        imgG = img.getGraphics();
        imgSize = sz;
        return true;
    }
    
    public void update(Graphics g) {
        if (img == null)
            g.clearRect(0, 0, getSize().width, getSize().height);
        paintComponent(g);
    }
    
    protected void paintComponent(Graphics g) {
        super.paintComponent(g);
        if (model != null) {
            newImgBuf(); // refresh the image buffer if necessary
            // compute the real-to-screen ratio, this variable differs
            // from model.real2Screen by zoomScale
            Dimension sz = getSize();
            double real2Screen0 = model.getScaleFromSpan(realSpan, sz.width, sz.height);
            model.setMatrix(viewMatrix, real2Screen0 * zoomScale,
                            sz.width/2, sz.height/2);
            imgG.setColor(Color.BLACK);
            imgG.fillRect(0, 0, sz.width, sz.height);
            model.paint(imgG);
            g.drawImage(img, 0, 0, this);
        }
    }
    
    /** Refresh the coordinates
     *  x[][] is the wrapped coordinates
     *  `n' may be less than x.length */
    public void refresh(double x[][], int n,
                        boolean center, boolean adjScale) {
        if (model == null) {
            model = new XYZModel();
            adjScale = true;
        }
        model.updateXYZ(x, n, center);
        if ( adjScale ) realSpan = model.getSpan(x, n);
        repaint();
    }
    
    /** Event handling */
    public void mouseClicked(MouseEvent e) { }
    public void mousePressed(MouseEvent e) {
        mouseX = e.getX();
        mouseY = e.getY();
        // consume this event so that it will not be processed
        // in the default manner
        e.consume();
    }
    public void mouseReleased(MouseEvent e) { }
    public void mouseEntered(MouseEvent e) { }
    public void mouseExited(MouseEvent e) { }
    
    public void mouseDragged(MouseEvent e) {
        int x = e.getX();
        int y = e.getY();
        tmpMatrix.unit();
        tmpMatrix.xrot(360.0 * (mouseY - y) / getSize().height);
        tmpMatrix.yrot(360.0 * (x - mouseX) / getSize().width);
        viewMatrix.mult(tmpMatrix);
        repaint();
        mouseX = x;
        mouseY = y;
        e.consume();
    }
    
    public void mouseMoved(MouseEvent e) { }
    
    public void mouseWheelMoved(MouseWheelEvent e) {
        int notches = e.getWheelRotation();
        if ((zoomScale -= 0.05f*notches) < 0.09999f)
            zoomScale = 0.1f;
        repaint();
    }
}


/** A set of atoms with 3D coordinates */
class XYZModel {
    double realXYZ[][]; // 3D real coordinates [np][3]
    int screenXYZ[][];  // 3D screen coordinates in pixels [np][3]
    // only the first two dimensions are used in drawing
    int zOrder[]; // z-order, larger is closer to the viewer
    int np = -1;
    
    boolean transformed;
    // rotation/scaling/translation matrix for the conversion from realXYZ[] to screenXYZ[]
    Matrix3D mat = new Matrix3D();
    double real2Screen = 50.0; // real size --> screen size
    double ballSize = 0.5; // ball size (radius) in terms of the real coordinates
    // 0.5 for hard spheres
    
    Atom atomDefault = new Atom(0.1, 0.7, 0.1, 1.0, 0.5, 0.5, 1.0);
    Atom atoms[]; // for colors of atoms
    
    XYZModel() {}
    
    /** Set the color of a particular atom */
    void setAtom(int i, Atom atom) {
        if ( i >= 0 && i < atoms.length )
            atoms[i] = atom;
    }
    
    /** Refresh coordinates
     *  x[0..n-1][3] is a three-dimensional vector
     *  n can be less than x.length */
    void updateXYZ(double x[][], int n, boolean center) {
        if (n != np) {
            //System.out.printf("XYZModel.updateXYZ n %d --> %d, (%g, %g, %g) (%g, %g, %g)\n", np, n, x[0][0], x[0][1], x[0][2], x[n-1][0], x[n-1][1], x[n-1][2]);
            np = n;
            realXYZ = new double [np][3];
            screenXYZ = new int [np][3];
            zOrder = new int [np];
            atoms = new Atom [np];
            for (int i = 0; i < np; i++)
                atoms[i] = atomDefault;
        }
        
        for (int d = 0; d < 3; d++) {
            double xc = 0;
            if ( center ) {
                for (int i = 0; i < np; i++)
                    xc += x[i][d];
                xc /= np;
            }
            for (int i = 0; i < np; i++)
                realXYZ[i][d] = x[i][d] - xc;
        }
        
        transformed = false;
    }
    
    /** Set the view matrix
     *  `s' is the scaling factor of translating real coordinates
     *    to the screen coordinates
     *  (x0, y0) the screen coordinates of the center */
    void setMatrix(Matrix3D viewMat, double s, double x0, double y0) {
        mat.unit();
        mat.mult(viewMat);
        mat.scale(s, s, s);
        real2Screen = s;
        mat.translate(x0, y0, 0);
        transformed = false;
    }
    
    /** Get the span of the model
     *  `n' may be less than x.length */
    double getSpan(double x[][], int n) {
        int dim = x[0].length;
        double realSpan = 0, del, fw, fh;
        
        for (int d = 0; d < dim; d++) {
            double xmin = 1e30, xmax = -1e30;
            for (int i = 0; i < n; i++)
                if ( x[i][d] < xmin ) xmin = x[i][d];
                else if ( x[i][d] > xmax ) xmax = x[i][d];
            if ( (del = xmax - xmin) > realSpan )
                realSpan = del;
        }
        return realSpan;
    }
    
    /** Translate a given span, return the real-to-screen ratio
     *  `w' and `h' are the width and height of the screen in pixels */
    double getScaleFromSpan(double realSpan, int w, int h) {
        realSpan += ballSize * 2; // add two radii
        double fw = w / realSpan;
        double fh = h / realSpan;
        double facShrink = 0.9; // shrink a bit for the margin
        return (fw < fh ? fw : fh) * facShrink;
    }
    
    /** Compute the Z-order */
    void getZOrder() {
        // transform the coordinates
        if ( !transformed ) {
            mat.transform(realXYZ, screenXYZ, np);
            transformed = true;
        }
        
        // bubble sort z-order
        // zOrder[0] is the fartherest from the viewer
        // zOrder[np - 1] is the nearest to the viewer
        for (int i = 0; i < np; i++)
            zOrder[i] = i;
        for (int i = 0; i < np; i++) {
            // find the particle with the smallest z
            int jm = i, k;
            for (int j = i + 1; j < np; j++)
                if (screenXYZ[zOrder[j]][2] < screenXYZ[zOrder[jm]][2])
                    jm = j;
            if (jm != i) {
                k = zOrder[i];
                zOrder[i] = zOrder[jm];
                zOrder[jm] = k;
            }
        }
    }
    
    /** Draw atom */
    void drawAtom(Graphics g, int iz) {
        int i = zOrder[iz];
        Atom atom = atoms[i];
        if (atom == null) return;
        double zMin = screenXYZ[zOrder[0]][2]; // fartherest from the viewer
        double zMax = screenXYZ[zOrder[np - 1]][2]; // closest to the viewer
        int greyScale = Atom.nBalls - 1;
        if (zMin != zMax) // zMin == zMax means the two-dimensional case
            greyScale = (int) (Atom.nBalls * (screenXYZ[i][2] - zMin) / (zMax - zMin) - 1e-6);
        // the atom closest to the viewer has a greyScale of Atom.nBalls - 1
        // the atom fartherest from the viewer has a greyScale of 0
        double radius = ballSize * atom.relRadius * real2Screen;
        atom.paint(g, screenXYZ[i][0], screenXYZ[i][1], greyScale, radius);
    }
    
    /** Paint this model to the graphics context `g' */
    void paint(Graphics g) {
        if (realXYZ == null || np <= 0) return;
        getZOrder();
        for (int iz = 0; iz < np; iz++)
            drawAtom(g, iz);
    }
}




/** Atom: a ball */
class Atom {
    private final static int R = 120;
    private final static int hx = 45; // (hx, hy) is the offset of the spot light from the center
    private final static int hy = 45;
    private static int maxr; // maximal distance from the spotlight
    public final static int nBalls = 16; // shades of grey
    // spotlight brightness (0.0, 1.0)
    // 1.0 means the spotlight is pure white
    // 0.0 means the spotlight is the same as the solid color
    double spotlightAmp = .4;
    // color damping along the distance from the spotlight
    double rContrast = 0.7;
    // z-depth contrast (0.0, inf)
    // inf means the fartherest atom is completely dark
    // 0.0 means the fartherest atom is the same as the
    //     nearest atom
    double zContrast = 2.0;
    
    private static byte data[];
    static { // data[] is a bitmap image of the ball of radius R
        data = new byte[R * 2 * R * 2];
        for (int Y = -R; Y < R; Y++) {
            int x0 = (int) (Math.sqrt(R * R - Y * Y) + 0.5);
            for (int X = -x0; X < x0; X++) {
                // sqrt(x^2 + y^2) gives distance from the spot light
                int x = X + hx, y = Y + hy;
                int r = (int) (Math.sqrt(x * x + y * y) + 0.5);
                // set the maximal intensity to the maximal distance
                // (in pixels) from the spot light
                if (r > maxr) maxr = r;
                    data[(Y + R) * (R * 2) + (X + R)] = (r <= 0) ? 1 : (byte) r;
                    }
        }
    }
    
    // the following variables are atom dependent
    private int Rl = 100, Gl = 100, Bl = 100;
    private Image balls[]; // 0..nBalls-1, at different z distances
    
    /** Constructor */
    Atom(int r, int g, int b) {
        setRGB(r, g, b);
    }
    
    /** colors that range from 0 to 1 */
    Atom(double r, double g, double b) {
        setRGB(r, g, b);
    }
    
    double relRadius = 1; // only used to save the radius
    
    Atom(double r, double g, double b, double rad) {
        this(r, g, b);
        relRadius = rad;
    }
    
    Atom(double r, double g, double b, double rad,
         double rcontrast, double spotlight, double zcontrast) {
        relRadius = rad;
        rContrast = rcontrast;
        spotlightAmp = spotlight;
        zContrast = zcontrast;
        setRGB(r, g, b);
    }
    
    /** Set color */
    void setRGB(int r, int g, int b) {
        Rl = r;
        Gl = g;
        Bl = b;
        makeBalls();
    }
    
    void setRGB(double r, double g, double b) {
        Rl = (int) (256*r - 1e-6);
        Gl = (int) (256*g - 1e-6);
        Bl = (int) (256*b - 1e-6);
        makeBalls();
    }
    
    /** Linearly interpolate colors */
    private int blend(double fg, double bg, double fgfactor) {
        return (int) (bg + (fg - bg) * fgfactor);
    }
    
    // need a component instance to call createImage()
    private static Component component = new Applet();
    
    /** Prepare ball images with different sizes */
    private void makeBalls() {
        balls = new Image[nBalls];
        byte red[]   = new byte[256];
        byte green[] = new byte[256];
        byte blue[]  = new byte[256];
        for (int id = 0; id < nBalls; id++) {
            // smaller `b' means closer to black
            // if id == 0 (fartherest from the viewer)
            //        b = 1/(1 + zContrast)
            //        the outer blend() gives a color close to bgGrey
            // if id == nBalls - 1 (closest to the viewer),
            //        b = 1, the outer blend() gives the color of
            //        the inner blend()
            double b = (zContrast*id/(nBalls - 1) + 1) / (zContrast + 1);
            for (int i = maxr; i >= 0; --i) {
                // closeness to the spotlight
                double q = 1 - 1. * i / maxr;
                // dampness of the color along the radius
                double p = 1 - rContrast * i / maxr;
                // contrast of the spotlight
                // if i == 0 (closest to the spotlight),
                //        d = 1.0 - spotlightAmp, the inner
                //        blend() gives a color close to 255
                //        (if spotlightAmp == 1).
                // if i == maxr (fartherest from the spotlight),
                //        d = 1.0, the inner blend() gives
                //        the foreground color, i.e., Rl, Gl, Bl
                // Thus, the inner blend() depends on the distance
                // from the spotlight, i == 0 means to be closest
                // to the spotlight
                double d = 1 - q * spotlightAmp;
                red[i]   = (byte) blend(blend(Rl*p, 255, d), 0, b);
                green[i] = (byte) blend(blend(Gl*p, 255, d), 0, b);
                blue[i]  = (byte) blend(blend(Bl*p, 255, d), 0, b);
            }
            // 256 color model
            IndexColorModel model = new IndexColorModel(
                                                        8, maxr + 1, red, green, blue, 0);
            balls[id] = component.createImage(
                                              new MemoryImageSource(R * 2, R * 2, model, data, 0, R * 2) );
        }
    }
    
    /** Draw a ball at screen coordinate (x, y) with a ball index `id'
     *  (0, 0) represents the top-left corner
     *  `x', `y', `radius' are given in pixels
     *  the ball index (gray code) `id' can be 0 to 15 */
    void paint(Graphics gc, int x, int y, int id, double radius) {
        if (balls == null) makeBalls();
        Image img = balls[id]; // id = [0..15]
        
        int size = (int) (radius * 2 + .5);
        gc.drawImage(img, x - size/2, y - size/2, size, size, null);
        //System.out.println("" + x + " " + y + " " + id + " " + radius);
    }
}

class Matrix3D {
    double xx, xy, xz, xo;
    double yx, yy, yz, yo;
    double zx, zy, zz, zo;
    static final double pi = 3.14159265;
    
    /** Create a new unit matrix */
    Matrix3D() {
        xx = 1.0;
        yy = 1.0;
        zz = 1.0;
    }
    
    /** Scale along each axis independently */
    void scale(double xf, double yf, double zf) {
        xx *= xf; xy *= xf; xz *= xf; xo *= xf;
        yx *= yf; yy *= yf; yz *= yf; yo *= yf;
        zx *= zf; zy *= zf; zz *= zf; zo *= zf;
    }
    
    /** Translate the origin */
    void translate(double x, double y, double z) {
        xo += x;
        yo += y;
        zo += z;
    }
    
    /** Rotate theta degrees around the y axis */
    void yrot(double theta) {
        theta *= (pi / 180);
        double ct = Math.cos(theta);
        double st = Math.sin(theta);
        
        double Nxx = (xx * ct + zx * st);
        double Nxy = (xy * ct + zy * st);
        double Nxz = (xz * ct + zz * st);
        double Nxo = (xo * ct + zo * st);
        
        double Nzx = (zx * ct - xx * st);
        double Nzy = (zy * ct - xy * st);
        double Nzz = (zz * ct - xz * st);
        double Nzo = (zo * ct - xo * st);
        
        xo = Nxo; xx = Nxx; xy = Nxy; xz = Nxz;
        zo = Nzo; zx = Nzx; zy = Nzy; zz = Nzz;
    }
    
    /** Rotate theta degrees about the x axis */
    void xrot(double theta) {
        theta *= (pi / 180);
        double ct = Math.cos(theta);
        double st = Math.sin(theta);
        
        double Nyx = (yx * ct + zx * st);
        double Nyy = (yy * ct + zy * st);
        double Nyz = (yz * ct + zz * st);
        double Nyo = (yo * ct + zo * st);
        
        double Nzx = (zx * ct - yx * st);
        double Nzy = (zy * ct - yy * st);
        double Nzz = (zz * ct - yz * st);
        double Nzo = (zo * ct - yo * st);
        
        yo = Nyo; yx = Nyx; yy = Nyy; yz = Nyz;
        zo = Nzo; zx = Nzx; zy = Nzy; zz = Nzz;
    }
    
    /** Rotate theta degrees about the z axis */
    void zrot(double theta) {
        theta *= pi / 180;
        double ct = Math.cos(theta);
        double st = Math.sin(theta);
        
        double Nyx = (yx * ct + xx * st);
        double Nyy = (yy * ct + xy * st);
        double Nyz = (yz * ct + xz * st);
        double Nyo = (yo * ct + xo * st);
        
        double Nxx = (xx * ct - yx * st);
        double Nxy = (xy * ct - yy * st);
        double Nxz = (xz * ct - yz * st);
        double Nxo = (xo * ct - yo * st);
        
        yo = Nyo; yx = Nyx; yy = Nyy; yz = Nyz;
        xo = Nxo; xx = Nxx; xy = Nxy; xz = Nxz;
    }
    
    /** Multiply this matrix by a second: M = M*R */
    void mult(Matrix3D rhs) {
        double lxx = xx * rhs.xx + yx * rhs.xy + zx * rhs.xz;
        double lxy = xy * rhs.xx + yy * rhs.xy + zy * rhs.xz;
        double lxz = xz * rhs.xx + yz * rhs.xy + zz * rhs.xz;
        double lxo = xo * rhs.xx + yo * rhs.xy + zo * rhs.xz + rhs.xo;
        
        double lyx = xx * rhs.yx + yx * rhs.yy + zx * rhs.yz;
        double lyy = xy * rhs.yx + yy * rhs.yy + zy * rhs.yz;
        double lyz = xz * rhs.yx + yz * rhs.yy + zz * rhs.yz;
        double lyo = xo * rhs.yx + yo * rhs.yy + zo * rhs.yz + rhs.yo;
        
        double lzx = xx * rhs.zx + yx * rhs.zy + zx * rhs.zz;
        double lzy = xy * rhs.zx + yy * rhs.zy + zy * rhs.zz;
        double lzz = xz * rhs.zx + yz * rhs.zy + zz * rhs.zz;
        double lzo = xo * rhs.zx + yo * rhs.zy + zo * rhs.zz + rhs.zo;
        
        xx = lxx; xy = lxy; xz = lxz; xo = lxo;
        yx = lyx; yy = lyy; yz = lyz; yo = lyo;
        zx = lzx; zy = lzy; zz = lzz; zo = lzo;
    }
    
    /** Recover the unit matrix */
    void unit() {
        xo = 0; xx = 1; xy = 0; xz = 0;
        yo = 0; yx = 0; yy = 1; yz = 0;
        zo = 0; zx = 0; zy = 0; zz = 1;
    }
    
    /** Transform np points from v into tv.
     *  v contains the input coordinates in floating point.
     *  Three successive entries in the array constitute a point.
     *  tv ends up holding the transformed points as integers;
     *  three successive entries per point */
    void transform(double v[][], int tv[][], int np) {
        // np can be different from v.length
        for (int i = 0; i < np; i++) {
            double x = v[i][0], y = v[i][1], z = v[i][2];
            tv[i][0] = (int) (xx * x + xy * y + xz * z + xo);
            tv[i][1] = (int) (yx * x + yy * y + yz * z + yo);
            tv[i][2] = (int) (zx * x + zy * y + zz * z + zo);
        }
    }
}



