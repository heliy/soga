package util;

public class StringArray {
    String[] array;
    
	public StringArray(String[] genos) {
		array=genos;
	}

	public int[] countAllChars(String[] letters) {
		// 统计每个字母出现的次数
     	int	len = letters.length;
		int[] counts = new int[len];
		for(String s: array){
			int i;
			System.out.println("+ " + s);
			for(i = 0;i < len;i++){
				if(letters[i] == s){
					counts[i]++;
					System.out.println("-" + letters[i] + " : "+ counts[i]);
					break;
				}
			}
		}
		return counts;
	}

}
