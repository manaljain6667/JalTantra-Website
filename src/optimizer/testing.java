package optimizer;
import java.util.*;
import java.text.*;
public class testing {
    public static void main(String[] args) {
        // Create a new file
        parseDate("00:05:00",16);
    }
    public static String parseDate(String dateString,int multiplier) {
        SimpleDateFormat dateFormat = new SimpleDateFormat("HH:mm:ss");
        String updatedTime="";
        try {
            System.out.println(dateString);
            Date date = dateFormat.parse(dateString);

            int doubled=date.getMinutes() * multiplier;
            System.out.println("minutes : "+ doubled);

            int hours = doubled / 60;
            int remainingMinutes = doubled % 60;
            int seconds = 0;

            updatedTime = String.format("%02d:%02d:%02d", hours, remainingMinutes, seconds);
            System.out.println("Time: " + updatedTime);


        } catch (ParseException e) {
            System.out.println("Error parsing date string: " + e.getMessage());
        }

        return updatedTime;
    }
}
/***
 File file = new File("sample.txt");
 if (!file.exists()) {
 try {
 file.createNewFile();
 // Write to a file
 String content = "1";
 try (FileWriter writer = new FileWriter(file)) {
 writer.write(content);
 System.out.println("File written successfully.");
 } catch (IOException e) {
 System.out.println("Error writing to file: " + e.getMessage());
 }
 System.out.println("File created successfully.");
 } catch (IOException e) {
 System.out.println("Error creating file: " + e.getMessage());
 }
 } else {
 // Read from a file
 try (BufferedReader reader = new BufferedReader(new FileReader(file))) {
 String line = reader.readLine();
 int attempt=Integer.parseInt(line);
 System.out.println("File content: " + attempt);
 } catch (IOException e) {
 System.out.println("Error reading from file: " + e.getMessage());
 }
 System.out.println("File already exists.");
 }


 public String parseDate(String dateString,int multiplier) {
 	SimpleDateFormat dateFormat = new SimpleDateFormat("HH:mm:ss");
 	String updatedTime="";
 	try {
 		System.out.println(dateString);
 		Date date = dateFormat.parse(dateString);

 		long millis = date.getMinutes();
 		System.out.println(millis + "MInutes qefsdfd");
 		long minutes = millis / 1000 / 60;

 		System.out.println("milliseconds "+millis+" minutes: "+minutes);

 		long doubledMinutes = 2 * minutes;
 		long doubledMillis = doubledMinutes * 60 * 1000;

 		Date doubledDate = new Date(doubledMillis);
 		System.out.println("Doubled time: " + dateFormat.format(doubledDate));
 		updatedTime=dateFormat.format(doubledDate);
 		System.out.println("Doubled time: " + updatedTime);

 	} catch (ParseException e) {
 		System.out.println("Error parsing date string: " + e.getMessage());
 	}

 	return updatedTime;
 }


 *
 */

