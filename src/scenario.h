class MyScenario {
  using World = MecaCell::BasicWorld<MyCell, MecaCell::Euler>;

 private:
  World w;

 public:
  World &getWorld() { return w; }

  void init(int argc, char** argv) {
    w.addCell(new Cell(MecaCell::Vec::zero()));
  }

  void loop(){
    w.update();
  }
};
