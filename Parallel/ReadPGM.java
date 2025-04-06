import java.io.*;
import java.util.Scanner;

public class ReadPGM {
    public static void read(double[][] X, String filename, int N) throws Exception {
        try (Scanner scanner = new Scanner(new File(filename))) {
            // Read PGM header
            String magicNumber = scanner.nextLine();
            if (!magicNumber.equals("P2") && !magicNumber.equals("P5")) {
                throw new IOException("Invalid PGM file format. Magic number must be P2 or P5");
            }

            // Skip comments
            String line = scanner.nextLine();
            while (line.startsWith("#")) {
                line = scanner.nextLine();
            }

            // Read dimensions
            String[] dimensions = line.split("\\s+");
            int width = Integer.parseInt(dimensions[0]);
            int height = Integer.parseInt(dimensions[1]);

            if (width != N || height != N) {
                throw new IOException("Image dimensions do not match expected size: " + N);
            }

            // Read max value
            int maxVal = scanner.nextInt();

            // Read pixel values
            for (int i = 0; i < N; i++) {
                for (int j = 0; j < N; j++) {
                    if (scanner.hasNextInt()) {
                        X[i][j] = (double) scanner.nextInt() / maxVal;
                    }
                }
            }
        } catch (IOException e) {
            System.err.println("Error reading PGM file: " + e.getMessage());
            throw e;
        }
    }
} 