
import java.util.Arrays ;

import java.awt.* ;
import javax.swing.* ;

public class Display2dFT extends JPanel {

    public static int CELL_SIZE = 1 ;

    int n ;

    double [] [] ftRe, ftIm ;

    Display2dFT(double [] [] ftRe, double [] [] ftIm, int n, String title) {

        setPreferredSize(new Dimension(CELL_SIZE * n, CELL_SIZE * n)) ;

        JFrame frame = new JFrame(title);
        frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        frame.setContentPane(this);
        frame.pack();
        frame.setVisible(true);

        this.ftRe = ftRe ;
        this.ftIm = ftIm ;

        this.n = n ;

        repaint() ;
    }

    public void paintComponent(Graphics g) {

        double [] mods = new double [n * n] ;
        for(int i = 0 ; i < n ; i++) {
            for(int j = 0 ; j < n ; j++) {
                mods [n * i + j] = mod(ftRe [i] [j], ftIm [i] [j]) ;
            }
        }
        Intensity intensity = new Intensity(mods) ;

        Hue hue = new Hue() ;

        for(int i = 0 ; i < n ; i++) {
            for(int j = 0 ; j < n ; j++) {

                // complex rendered as a maximally saturated colour,
                // brightness dependent on modulus, and hue on
                // red-green-blue color wheel determined by argument.

                int iSigned = i < n/2 ? i : i - n ;
                int jSigned = j < n/2 ? j : j - n ;
                double re = ftRe [n/2 + iSigned] [n/2 - 1 - jSigned] ;

                double im = ftIm [n/2 + iSigned] [n/2 - 1 - jSigned] ;

                hue.set(re, im) ;
                double mod = mod(re, im) ;

                double log = Math.log(mod) ;
                double brightness = intensity.mul * log + intensity.con ;
                if(brightness < 0) {
                    brightness = 0 ;
                }

                float red = (float) (brightness * hue.red) ;
                float green = (float) (brightness * hue.green) ;
                float blue = (float) (brightness * hue.blue) ;

                Color c = new Color(red, green, blue) ;
                g.setColor(c) ;
                g.fillRect(CELL_SIZE * i, CELL_SIZE * j,
                        CELL_SIZE, CELL_SIZE) ;
            }
        }
    }

    static class Intensity {

        // Scaling for colour intensity based on modulus of complex number.

        // An underlying logarithmic scale is assumed.

        public static final double MEDIAN_MOD_BRIGHTNESS = 0.3 ;

        double mul, con ;  // "m" and "c" parameters in linear plot

        Intensity(double [] mods) {

            int n = mods.length ;

            Arrays.sort(mods, 0, n) ;

            // Maximum and median moduli
            double modMax = mods [n - 1] ;
            double modMed = 0.5 * (mods [n / 2 - 1] + mods [n / 2]) ;

            double logModMax = Math.log(modMax) ;
            double logModMed = Math.log(modMed) ;

            // m and c in linear equation that maps max mod to 1.0 and
            // median mod to MEDIAN_MOD_BRIGHTNESS:

            mul = (1 - MEDIAN_MOD_BRIGHTNESS) / (logModMax - logModMed) ;

            con = 1 - mul * logModMax ;
        }
    }

    static class Hue {

        // Associate a "hue" with argument of a complex number, according
        // to some prescription (red is centred on positive real axis,
        // green at angle 2pi/3 in complex plane, blue at 4pi/3).

        // red + green + blue = 1.0

        double red, green, blue ;

        void set(double re, double im) {

            double arg = arg(re, im) ;

            if(arg < 0) {
                arg += 2 * Math.PI ;  // makes logic below simpler
            }
            double argNorm = 3 * arg / (2 * Math.PI) ;
            // in range 0 to 3

            if(argNorm < 1) {
                green = argNorm ;
                red = 1 - red ;
                blue = 0 ;
            }
            else if(argNorm < 2) {
                blue = argNorm - 1 ;
                green = 1 - blue ;
                red = 0 ;
            }
            else {
                red = argNorm - 2 ;
                blue = 1 - red ;
                green = 0 ;
            }
        }
    }

    // Complex support

    static double arg(double re, double im) {
        double arg ;
        if(re == 0) {
            if(im >= 0) {
                arg = Math.PI / 2 ;
            }
            else {
                arg = - Math.PI / 2 ;
            }
        }
        else {
            arg = Math.atan(im / re) ;
            if(re < 0) {
                if(arg > 0) {
                    arg -= Math.PI ;
                }
                else {
                    arg += Math.PI ;
                }
            }
        }
        return arg ;
    }

    static double mod(double re, double im) {

        return Math.sqrt(re * re + im * im) ;
    }

}
