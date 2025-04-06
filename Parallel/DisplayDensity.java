import java.awt.*;
import javax.swing.*;

public class DisplayDensity extends JPanel {
    
    public static int CELL_SIZE = 1;
    
    int n;
    double[][] data;
    
    public DisplayDensity(double[][] data, int n, String title) {
        this.data = data;
        this.n = n;
        
        setPreferredSize(new Dimension(CELL_SIZE * n, CELL_SIZE * n));
        
        JFrame frame = new JFrame(title);
        frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        frame.setContentPane(this);
        frame.pack();
        frame.setVisible(true);
        
        repaint();
    }
    
    @Override
    public void paintComponent(Graphics g) {
        super.paintComponent(g);
        
        // Find min and max values for normalization
        double min = Double.MAX_VALUE;
        double max = Double.MIN_VALUE;
        
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                min = Math.min(min, data[i][j]);
                max = Math.max(max, data[i][j]);
            }
        }
        
        // Normalize and draw
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                // Normalize the value to [0, 1]
                float normalizedValue;
                if (max > min) {
                    normalizedValue = (float) ((data[i][j] - min) / (max - min));
                } else {
                    normalizedValue = 0.5f;
                }
                
                // Convert to grayscale color
                int grayValue = (int) (normalizedValue * 255);
                Color color = new Color(grayValue, grayValue, grayValue);
                
                g.setColor(color);
                g.fillRect(i * CELL_SIZE, j * CELL_SIZE, CELL_SIZE, CELL_SIZE);
            }
        }
    }
} 