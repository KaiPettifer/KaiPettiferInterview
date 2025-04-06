import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;

public class MultiFTT {

    public static int N = 256;
    //public static int Nthreads = Runtime.getRuntime().availableProcessors();
    public static int Nthreads = 1;

    public static void main(String[] args) throws Exception {

        double[][] X = new double[N][N];
        ReadPGM.read(X, "wolf.pgm", N);

        //display the first image before any changes are made
        DisplayDensity display = new DisplayDensity(X, N, "Original Image");

        //start time to compute time taken
        long startDFT = System.currentTimeMillis();

        double[][] CRe = new double[N][N];
        double[][] CIm = new double[N][N];
        for (int k = 0; k < N; k++) {
            System.arraycopy(X[k], 0, CRe[k], 0, N);
        }

         // Fourier transform and end time once complete and prints difference in time
        fft2d(CRe, CIm, 1);
        long endDFT = System.currentTimeMillis();
        System.out.println("DFT computation time: " + (endDFT - startDFT) + " ms");
        Display2dFT display2 = new Display2dFT(CRe, CIm, N, "Discrete FT");

        // start another timer for inverse
        long startInverse = System.currentTimeMillis();
        double[][] reconRe = new double[N][N];
        double[][] reconIm = new double[N][N];
        for (int k = 0; k < N; k++) {
            System.arraycopy(CRe[k], 0, reconRe[k], 0, N);
            System.arraycopy(CIm[k], 0, reconIm[k], 0, N);
        }

        //inverse to reverse image and end timer to print time taken and inverse image
        fft2d(reconRe, reconIm, -1);
        long endInverse = System.currentTimeMillis();
        DisplayDensity display3 = new DisplayDensity(reconRe, N, "Reconstructed Image");
        System.out.println("Inverse computation time: " + (endInverse - startInverse) + " ms");
    }

   //define method introducing threads
    static void fft2d(double[][] re, double[][] im, int isgn) {
        ExecutorService rowExecutor = Executors.newFixedThreadPool(Nthreads);

        // loops through each row handing each row to a different threads
        for (int i = 0; i < N; i++) {
            final int row = i;
            rowExecutor.execute(() -> fft1d(re[row], im[row], isgn));
        }

        rowExecutor.shutdown();
        try {
            rowExecutor.awaitTermination(Long.MAX_VALUE, TimeUnit.NANOSECONDS);
        } catch (InterruptedException e) {
            e.printStackTrace();
        }

        ExecutorService colExecutor = Executors.newFixedThreadPool(Nthreads);

        // Transpose
        double[][] transposedRe = transpose(re);
        double[][] transposedIm = transpose(im);

        // loops each row again for each transposed row to perform fft1d
        for (int i = 0; i < N; i++) {
            final int row = i;
            colExecutor.execute(() -> fft1d(transposedRe[row], transposedIm[row], isgn));
        }

        colExecutor.shutdown();
        try {
            colExecutor.awaitTermination(Long.MAX_VALUE, TimeUnit.NANOSECONDS);
        } catch (InterruptedException e) {
            e.printStackTrace();
        }

        // Transpose back
        transpose(transposedRe);
        transpose(transposedIm);


        for (int k = 0; k < N; k++) {
            System.arraycopy(transposedRe[k], 0, re[k], 0, N);
            System.arraycopy(transposedIm[k], 0, im[k], 0, N);
        }
    }

    //define transpose 
    static double[][] transpose(double[][] a) {
        int rows = a.length;
        int cols = a[0].length;

        double[][] transposed = new double[cols][rows];
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                transposed[j][i] = a[i][j];
            }
        }

        return transposed;
    }

    //define fft1d which ises bit reverse to compute each row and not 2d
    static void fft1d(double[] re, double[] im, int isgn) {
        int N = re.length;

        // Bit-reverse the input
        int[] bitrev = bitrev(N);
        for (int k = 0; k < N; k++) {
            int j = bitrev[k];
            if (j > k) {
                double temp = re[k];
                re[k] = re[j];
                re[j] = temp;
                temp = im[k];
                im[k] = im[j];
                im[j] = temp;
            }
        }

        // Compute the FFT
        for (int n = 2; n <= N; n *= 2) {
            for (int i = 0; i < N; i += n) {
                for (int j = 0; j < n / 2; j++) {
                    int k1 = i + j;
                    int k2 = i + j + n / 2;
                    double w_re = Math.cos(2 * Math.PI * j / n);
                    double w_im = isgn * Math.sin(2 * Math.PI * j / n);
                    double t_re = w_re * re[k2] - w_im * im[k2];
                    double t_im = w_re * im[k2] + w_im * re[k2];
                    re[k2] = re[k1] - t_re;
                    im[k2] = im[k1] - t_im;
                    re[k1] += t_re;
                    im[k1] += t_im;
                }
            }
        }

        // Scale the output for the inverse FFT
        if (isgn < 0) {
            for (int i = 0; i < N; i++) {
                re[i] /= N;
                im[i] /= N;
            }
        }
    }

    static int[] bitrev(int N) {
        int[] bitrev = new int[N];

        for (int k = 0; k < N; k++) {
            bitrev[k] = Integer.reverse(k) >>> (32 - Integer.numberOfTrailingZeros(N));
        }

        return bitrev;
    }
}

