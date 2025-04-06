import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;

public class FutureFTT{

    public static int N = 256;
    //public static final int Nthreads = Runtime.getRuntime().availableProcessors();
    public static final int Nthreads = 1;

    public static void main(String[] args) throws Exception {
        double[][] X = new double[N][N];
        ReadPGM.read(X, "wolf.pgm", N);

        DisplayDensity display = new DisplayDensity(X, N, "Original Image");
        //get time to measure time difference between start and end
        long startDFT = System.currentTimeMillis();


        double[][] CRe = new double[N][N], CIm = new double[N][N];
        for (int k = 0; k < N; k++) {
            System.arraycopy(X[k], 0, CRe[k], 0, N);
        }

        //run transformaton and measure time to compare with first time
        fft2d(CRe, CIm, 1);
        long endDFT = System.currentTimeMillis();
        System.out.println("DFT computation time: " + (endDFT - startDFT) + " ms");
        Display2dFT display2 = new Display2dFT(CRe, CIm, N, "Discrete FT");

        // create array for in-place inverse FFT, and copy FT to it
        long startInverse = System.currentTimeMillis();
        double[][] reconRe = new double[N][N], reconIm = new double[N][N];
        for (int k = 0; k < N; k++) {
            System.arraycopy(CRe[k], 0, reconRe[k], 0, N);
            System.arraycopy(CIm[k], 0, reconIm[k], 0, N);
        }


        fft2d(reconRe, reconIm, -1); 
        long endInverse = System.currentTimeMillis();
        DisplayDensity display3 = new DisplayDensity(reconRe, N, "Reconstructed Image");
        System.out.println("Inverse computation time: " + (endInverse - startInverse) + " ms");
    }

    static void fft2d(double[][] re, double[][] im, int isgn) {
        //creates thread pools using executor and future to assign rows to each thread
        ExecutorService executor = Executors.newFixedThreadPool(Nthreads);
        Future<?>[] futures = new Future<?>[N];

        // Perform FFT on all rows of the input arrays
        for (int i = 0; i < N; i++) {
            final int row = i;
        //assigns each row to future array using executor
            futures[i] = executor.submit(() -> fft1d(re[row], im[row], isgn));
        }

        // Wait for all row FFTs to complete
        for (Future<?> future : futures) {
            try {
                future.get();
            } catch (Exception e) {
                e.printStackTrace();
            }
        }

        // Transpose
        transpose(re);
        transpose(im);

        // Perform FFT on all rows of the transposed arrays
        for (int i = 0; i < N; i++) {
            final int row = i;
            futures[i] = executor.submit(() -> fft1d(re[row], im[row], isgn));
        }

        // Wait for all transposed row FFTs to complete
        for (Future<?> future : futures) {
            try {
                future.get();
            } catch (Exception e) {
                e.printStackTrace();
            }
        }

        // Transpose the arrays back to their original orientation
        transpose(re);
        transpose(im);

        executor.shutdown();
    }

    static void transpose(double[][] a) {
        int N = a.length;

        for (int i = 0; i < N; i++) {
            for (int j = 0; j < i; j++) {
                // Swap the elements at indices (i,j) and (j,i)
                double temp = a[i][j];
                a[i][j] = a[j][i];
                a[j][i] = temp;
            }
        }
    }

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