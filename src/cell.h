class MyCell : public MecaCell::ConnectableCell<MyCell> {

  double getAdhesionWith(const MyCell *c) { return 0.9; }

  MyCell* updateBehavior(double deltaTime) {
    return nullptr;
  }
};
