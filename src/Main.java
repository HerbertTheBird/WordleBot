import java.io.*;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.*;
import java.util.stream.IntStream;

public class Main {
	static String[] guesses;
	static String[] answers;
	final static int L = 5;
	static int[] pow3 = new int[L];
	final static int sz = (int)Math.pow(3, L);
	static int[][] table;
	static double[] LOG2;
	static final double INV_LN2 = 1.0 / Math.log(2);

	static final ThreadLocal<int[]> TL_BUCKET  = ThreadLocal.withInitial(() -> new int[sz]);
	static final ThreadLocal<int[]> TL_TOUCHED = ThreadLocal.withInitial(() -> new int[sz]);
	
	static HashMap<StateKey, Integer> solutions = new HashMap();
	
	public static void main(String[] args) throws Exception{
		long start = System.currentTimeMillis();
		setup();
		System.out.println("setup takes " + (System.currentTimeMillis()-start) + " ms");
		int[] remaining = new int[answers.length];
		for(int i = 0; i < answers.length; i++) remaining[i] = i;
//		solve();
		Scanner in = new Scanner(System.in);
		System.out.println(guesses[solutions.get(new StateKey(remaining))]);
		int prevBest = 10183;
		while(remaining.length > 1) {
			int idx = -1;
			do {
				System.out.println("What word did you use?");
				String word = in.nextLine().toLowerCase();
				if(word.length() == 0) {
					idx = prevBest;
				}
				for(int i = 0; i < guesses.length; i++) {
					if(guesses[i].equals(word)) {
						idx = i;
						break;
					}
				}
			}while(idx == -1);
			int p = -1;
			do {
				System.out.println("What is the score? (2 = green, 1 = yellow, 0 = grey)");
				String pattern = in.nextLine();
				if(pattern.length() != 5) continue;
				p = 0;
				for(int i = 0; i < L; i++) {
					p += (pattern.charAt(i)-'0')*pow3[i];
					if(pattern.charAt(i) > '2' || pattern.charAt(i) < '0') {
						p = -1;
						break;
					}
				}
			}while(p == -1);
			int cnt = 0;
			for(int i = 0; i < remaining.length; i++) {
				if(table[idx][remaining[i]] == p) {
					cnt ++;
				}
			}
			int[] newrem = new int[cnt];
			cnt = 0;
			for(int i = 0; i < remaining.length; i++) {
				if(table[idx][remaining[i]] == p) {
					newrem[cnt++] = remaining[i];
				}
			}
			remaining = newrem;
			System.out.println(guesses[solutions.get(new StateKey(remaining))]);
		}
		
	}
	public static void solve() throws Exception{
		int[] remaining = new int[answers.length];
		for(int i = 0; i < answers.length; i++) remaining[i] = i;
		long start = System.currentTimeMillis();
		Pair[] a = go(remaining, 0.7, 32, 6); //solve - 30 seconds
		createHashMapFile();
		createAnswerFile();
		System.out.println("search takes " + (System.currentTimeMillis()-start) + " ms");
		for(int i = 0; i < a.length; i++) {
			System.out.println(guesses[a[i].s] + " " + a[i].d);
		}
	}
	public static void setup() throws Exception{
		pow3[0] = 1;
		for(int i = 1; i < L; i++)
			pow3[i] = 3*pow3[i-1];
		guesses = Files.readAllLines(Paths.get("guesses.txt")).toArray(new String[0]);
		answers = Files.readAllLines(Paths.get("answers.txt")).toArray(new String[0]);
		table = new int[guesses.length][answers.length];
		LOG2 = new double[answers.length + 1];
	    LOG2[0] = 0.0;
	    for (int k = 1; k <= answers.length; k++) LOG2[k] = Math.log(k) * INV_LN2;
		FileReader in = new FileReader("table.txt");
		for(int i = 0; i < table.length; i++) {
			for(int j = 0; j < table[0].length; j++) {
				table[i][j] = in.read();
				if(table[i][j] != score(guesses[i], answers[j])) {
					System.out.println("mismatch " + table[i][j] + " " + score(guesses[i], answers[j]) + " " + (char)table[i][j]);
				}
			}
		}
		in.close();
		BufferedReader in2 = new BufferedReader(new FileReader("solutionhm.txt"));
		String line;
		StringTokenizer st;
		while((line = in2.readLine()) != null) {
			st = new StringTokenizer(line);
			int n = Integer.parseInt(st.nextToken());
			int[] a = new int[n];
			for(int i = 0; i < n; i++) {
				a[i] = Integer.parseInt(st.nextToken());
			}
			int word = Integer.parseInt(st.nextToken());
			solutions.put(new StateKey(a), word);
		}
		in2.close();
	}
	public static void createHashMapFile() throws Exception{
		FileWriter writer = new FileWriter("solutionhm.txt");
		for(Map.Entry<StateKey, Integer> i:solutions.entrySet()) {
			writer.write(i.getKey().a.length + " ");
			for(int j = 0; j < i.getKey().a.length; j++) {
				writer.write(i.getKey().a[j] + " ");
			}
			writer.write(" " + i.getValue() + "\n");
		}
		writer.close();
	}
	public static void createAnswerFile() throws Exception{
		FileWriter writer = new FileWriter("solution.txt");
		for(int ans = 0; ans < answers.length; ans++) {
			int[] remaining = new int[answers.length];
			for(int j = 0; j < answers.length; j++) remaining[j] = j;
			while(true) {
				int s = solutions.get(new StateKey(remaining));
				System.out.print(guesses[s] + " ");
				int p = table[s][ans];
				writer.write(guesses[s]);
				int cnt = 0;
				for(int i = 0; i < remaining.length; i++) {
					if(table[s][remaining[i]] == p) {
						cnt ++;
					}
				}
				int[] newrem = new int[cnt];
				cnt = 0;
				for(int i = 0; i < remaining.length; i++) {
					if(table[s][remaining[i]] == p) {
						newrem[cnt++] = remaining[i];
					}
				}
				remaining = newrem;
				if(s == ans) break;
				else {
					writer.write(",");
				}
			}
			System.out.println();
			writer.write("\n");
		}
		writer.close();
	}
	public static double guessEntropy(int gIdx, int[] remaining) {
		int[] bucket = TL_BUCKET.get();
		int[] touched = TL_TOUCHED.get();
		int R = remaining.length;
		
        int tSz = 0;

        for (int k = 0; k < R; k++) {
            int a = remaining[k];
            int pat = table[gIdx][a];
            if (bucket[pat]++ == 0) touched[tSz++] = pat;
        }

        double sum = 0.0;
        for (int i = 0; i < tSz; i++) {
            int c = bucket[touched[i]];
            sum += c * LOG2[c];
            if(touched[i] == sz-1) {
            	sum -= 0.00001;
            }
        }
        double H = LOG2[R] - sum / R;
		for (int i = 0; i < tSz; i++) bucket[touched[i]] = 0;
		return H;
	}
	public static Pair[] go(int[] remaining, double sens, int K, int depth) {
		if(remaining.length == 1) {
			solutions.put(new StateKey(remaining), remaining[0]);
			return new Pair[]{new Pair(remaining[0], 1)}; //first few in guess are answers
		}
		final int R = remaining.length;


		PriorityQueue<Pair> pq = IntStream.range(0, guesses.length).parallel()
		    .mapToObj(gIdx -> new Pair(gIdx, guessEntropy(gIdx, remaining)))
		    .collect(
		        () -> new PriorityQueue<>(Comparator.comparingDouble(p -> p.d)),
		        (q, p) -> {
		            if (q.size() < K) q.offer(p);
		            else if (p.d > q.peek().d) {
		                q.poll();
		                q.offer(p);
		            }
		        },
		        (q1, q2) -> {
		            for (Pair p : q2) {
		                if (q1.size() < K) q1.offer(p);
		                else if (p.d > q1.peek().d) {
		                    q1.poll();
		                    q1.offer(p);
		                }
		            }
		        }
		    );

		double best = pq.stream().mapToDouble(p -> p.d).max().orElse(Double.NEGATIVE_INFINITY);
		double cutoff = best - sens;

		List<Pair> outList = new ArrayList<>();
		for (Pair p : pq) if (p.d >= cutoff) outList.add(p);
		outList.sort((a, b) -> Double.compare(b.d, a.d));
		if(remaining.length == answers.length) {
			outList = new ArrayList<Pair>();
			outList.add(new Pair(10183, 0));
		}
		Pair[] out = outList.toArray(new Pair[0]);

		if(depth == 0) {
			System.out.println("approximation");
			return out;
		}
		
		Arrays.parallelSetAll(out, t -> {
		    Pair p = out[t];

		    int[] cnt = new int[sz];
		    int[] touched = new int[sz];
		    int tSz = 0;
		    for (int k = 0; k < R; k++) {
		        int a = remaining[k];
		        int pat = table[p.s][a] & 0xFF;
		        if (cnt[pat]++ == 0) touched[tSz++] = pat;
		    }

		    int[] patToIdx = new int[sz];
		    for (int i = 0; i < tSz; i++) patToIdx[touched[i]] = i;

		    int[][] groups = new int[tSz][];
		    for (int i = 0; i < tSz; i++) groups[i] = new int[cnt[touched[i]]];

		    Arrays.fill(cnt, 0);

		    for (int k = 0; k < R; k++) {
		        int a = remaining[k];
		        int pat = table[p.s][a] & 0xFF;
		        int gi = patToIdx[pat];
		        groups[gi][ cnt[pat]++ ] = a;
		    }

		    double childExp = 0.0;
		    for (int i = 0; i < tSz; i++) {
		        int c = groups[i].length;
		        double prob = (double) c / R;
		        double d;
		        if(touched[i] == sz-1) {
		        	d = 1;
		        }
		        else {
		        	Pair childBest = go(groups[i], sens, K, depth - 1)[0];
		        	d = childBest.d+1;
		        }
		        childExp += prob * d;
		    }

		    return new Pair(p.s, childExp);
		});
		Arrays.sort(out);
		solutions.put(new StateKey(remaining), out[0].s);
		return out;
	}
	public static void createTable() throws Exception{
		FileWriter writer = new FileWriter("table.txt");
		for(int i = 0; i < guesses.length; i++) {
			for(int j = 0; j < answers.length; j++) {
				writer.write((char)score(guesses[i], answers[j]));
			}
		}
		writer.close();
	}
	public static int score(String guess, String answer) {
		int out = 0;
		int[] has = new int[26];
		for(int i = 0; i < L; i++)
			has[answer.charAt(i)-'a']++;
		for(int i = 0; i < L; i++) {
			if(guess.charAt(i) == answer.charAt(i)) {
				out += 2*pow3[i];
				has[guess.charAt(i)-'a']--;
			}
		}
		for(int i = 0; i < L; i++) {
			if((out / pow3[i])%3 != 2 && has[guess.charAt(i)-'a'] > 0) {
				out += pow3[i];
				has[guess.charAt(i)-'a']--;
			}
		}
		return out;
	}
}
class Pair implements Comparable<Pair>{
    public int s;
    public double d;

    public Pair(int s, double d) {
        this.s = s;
        this.d = d;
    }
    @Override
    public String toString() {
        return "(" + s + ", " + d + ")";
    }
    @Override
    public int compareTo(Pair other) {
    	int cmp = Double.compare(this.d, other.d);
        if (cmp != 0) return cmp;
        return Integer.compare(s, other.s);
    }
}

final class StateKey {
	  final int[] a;
	  final int hash;
	  StateKey(int[] rem) {
	    a = Arrays.copyOf(rem, rem.length);
	    Arrays.sort(a);
	    hash = Arrays.hashCode(a);
	  }
	  @Override public int hashCode() { return hash; }
	  @Override public boolean equals(Object o) {
	    return (o instanceof StateKey) && Arrays.equals(a, ((StateKey)o).a);
	  }
	}