struct BcData
{
    const double U1;
    const double U2;

    const double t_load;
    const double tg_load;

    //    static std::string inflow_data;

    const double wall;
    const double d_in;
    const double d_out;
    const double h1_in;
    const double h2_in;
    const double h_out;
    const double z_min;
    const double z_max;

    BcData(Pointer<Database> input_db)
        : U1(input_db->getDouble("U1")),
          U2(input_db->getDouble("U2")),
          t_load(input_db->getDouble("T_LOAD")),
          tg_load(input_db->getDouble("TG_LOAD")),
          wall(input_db->getDouble("WALL")),
          d_in(input_db->getDouble("D_IN")),
          d_out(input_db->getDouble("D_OUT")),
          h1_in(input_db->getDouble("H1_IN")),
          h2_in(input_db->getDouble("H2_IN")),
          h_out(input_db->getDouble("H_OUT")),
          z_min(input_db->getDouble("Z_MIN")),
          z_max(input_db->getDouble("Z_MAX"))
    {
    }
};
