import java.util.*;

final int D = 8;
final int DIAMETER = 24;
final int TS = 42;
final int ARROW = 150;
final double TA = 0.2;

String[] texts = {
  " 4-fold  D₃, A₃, B₂",
  " 6-fold  D₄, B₃, G₂, A₂",
  " 6-fold  A₅",
  " 6-fold  Lisi Triality  (12/2)",
  " 8-fold  D₅, B₄",
  " 8-fold  A₇",
  " 8-fold  B₈  (16/2)",
  "10-fold  D₆, B₅, A₄",
  "12-fold  D₇, B₆",
  "12-fold  E₆, F₄",
  "14-fold  D₈, B₇, A₆",
  "18-fold  E₇",
  "20-fold",
  "24-fold",
  "30-fold  E₈, H₄  Coxeter plane"
};

IntList crystals = new IntList();
ArrayList <Boolean> labels = new ArrayList <Boolean> ();
ArrayList <Double> scale = new ArrayList <Double> ();
ArrayList <Double [] [] > data = new ArrayList <Double[] [] > ();


color[] colors = {
  color(255, 255, 255),
  color(255, 165, 175),
  color(100, 255, 125),
  color(180, 175, 255),
  color(225, 100, 235)
};

final int[] [] ms = {
  {1, 12, 32, 60},
  {1, 27, 72},
  {3, 8, 24, 30},
  {1, 2, 3},
  {1, 8, 24},
  {1, 2, 4, 8},
  {1},
  {1, 5, 10, 20},
  {1, 3, 9, 12},
  {1, 8, 24},
  {1, 2, 3},
  {1, 3, 6},
  {1},
  {1},
  {1}
};

void setup() {
  size(1920, 120);
  frameRate(300);
  
  noStroke();
  ellipseMode(CENTER);
    
  file();
}

void file() {
  BufferedReader reader = createReader("/home/sam/Research/E8/QC2/data.tsv");
  String line = null;
  
  try {
    while ((line = reader.readLine()) != null) {
      String[] s = split(line, TAB);
      
      int     f     = int(s[0] );
      int     n     = int(s[1] );
      int     ff    = int(s[2] );
      boolean label = s[3].equals("0");
      double  rot   = Double.parseDouble(s[4] );
      double  zoom  = Double.parseDouble(s[5] );
      
      crystals.append(n);
      labels.add(label);
      scale.add(zoom);
      
      Double[] [] xy = new Double[2] [8];
      
      line = reader.readLine();
      s = split(line, TAB);
      
      for (int d = 0; d < D; d++)
        xy[0] [d] = Double.parseDouble(s[d] );
      
      line = reader.readLine();
      s = split(line, TAB);

      for (int d = 0; d < D; d++)
        xy[1] [d] = Double.parseDouble(s[d] );

      data.add(xy);
    }

    reader.close();
  } catch (IOException e) {
    e.printStackTrace();
  }
}

void draw() {
  if (frameCount >= data.size()) {
    exit();
    return; 
  }
  
  background(0, 0, 0);

  fill(255, 255, 255);
  
  final int     c = crystals.get(frameCount);
  final boolean l = labels  .get(frameCount);
    
  if (l) {
    textSize(1.25*TS);
    int TV = (int) (height/2 - TA*textAscent());
    
    textAlign(CENTER, CENTER);
    text(texts[c    ], width/4          , TV);
    
    textAlign(  LEFT, CENTER);
    
    textSize(1.5*TS);
    TV = (int) (height/2 - TA*textAscent());
    
    final int shift = c == 2 ? 1 : 0;
    final int len = ms[c].length;
    
    for (int i = 0; i < len; i++) {
      fill(colors[i + shift] );
      ellipse(width/2 + 120 + 210*i + (4 - len)*210/2, height/2, DIAMETER, DIAMETER); 
      text(ms[c] [i], width/2 + 120 + 210*i + 40 + (4 - len)*210/2, TV);
    }
  } else {
    textSize(1.25*TS);
    int TV = (int) (height/2 - TA*textAscent());
    
    textAlign( RIGHT, CENTER);
    text(texts[c    ], width/2 - ARROW/2, TV);
    textAlign(  LEFT, CENTER);
    text(texts[c + 1], width/2 + ARROW/2, TV);
    
    rect(width/2 - 0.7*ARROW/2, height/2 - 15, 0.7*ARROW - 15, 30);
    
    triangle(
      width/2 + 0.7*ARROW/2 + 12, height/2     ,
      width/2 + 0.7*ARROW/2 - 30, height/2 - 35,
      width/2 + 0.7*ARROW/2 - 30, height/2 + 35
    );
  }

  saveFrame("/home/sam/Research/E8/FT/Stills2/#####.png");
}
