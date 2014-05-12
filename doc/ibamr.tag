<?xml version='1.0' encoding='UTF-8' standalone='yes' ?>
<tagfile>
  <compound kind="class">
    <name>IBAMR::AdvDiffCenteredConvectiveOperator</name>
    <filename>class_i_b_a_m_r_1_1_adv_diff_centered_convective_operator.html</filename>
    <base>IBAMR::ConvectiveOperator</base>
    <member kind="function">
      <type></type>
      <name>AdvDiffCenteredConvectiveOperator</name>
      <anchorfile>class_i_b_a_m_r_1_1_adv_diff_centered_convective_operator.html</anchorfile>
      <anchor>a9a234f0435ac41caf320956e87f9a077</anchor>
      <arglist>(const std::string &amp;object_name, SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::CellVariable&lt; NDIM, double &gt; &gt; Q_var, SAMRAI::tbox::Pointer&lt; SAMRAI::tbox::Database &gt; input_db, ConvectiveDifferencingType difference_form, const std::vector&lt; SAMRAI::solv::RobinBcCoefStrategy&lt; NDIM &gt; * &gt; &amp;bc_coefs)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~AdvDiffCenteredConvectiveOperator</name>
      <anchorfile>class_i_b_a_m_r_1_1_adv_diff_centered_convective_operator.html</anchorfile>
      <anchor>aae1cbeffeaa0a1b43558f1b9e012c920</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>applyConvectiveOperator</name>
      <anchorfile>class_i_b_a_m_r_1_1_adv_diff_centered_convective_operator.html</anchorfile>
      <anchor>a9860b20c012d1b516418f5469fc8a88c</anchor>
      <arglist>(int Q_idx, int N_idx)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>initializeOperatorState</name>
      <anchorfile>class_i_b_a_m_r_1_1_adv_diff_centered_convective_operator.html</anchorfile>
      <anchor>a4348518ed962e1df23ffb55c0a3fe32c</anchor>
      <arglist>(const SAMRAI::solv::SAMRAIVectorReal&lt; NDIM, double &gt; &amp;in, const SAMRAI::solv::SAMRAIVectorReal&lt; NDIM, double &gt; &amp;out)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>deallocateOperatorState</name>
      <anchorfile>class_i_b_a_m_r_1_1_adv_diff_centered_convective_operator.html</anchorfile>
      <anchor>a7cfce9ccd0c7d1b5adb663f935e11617</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>ConvectiveOperator</name>
      <anchorfile>class_i_b_a_m_r_1_1_convective_operator.html</anchorfile>
      <anchor>a7ff1d2d0591334559ae2d428e8beda15</anchor>
      <arglist>(const std::string &amp;object_name, ConvectiveDifferencingType difference_form)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~ConvectiveOperator</name>
      <anchorfile>class_i_b_a_m_r_1_1_convective_operator.html</anchorfile>
      <anchor>a21b1ed229553121d8bc481fea89a9000</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setAdvectionVelocity</name>
      <anchorfile>class_i_b_a_m_r_1_1_convective_operator.html</anchorfile>
      <anchor>a92de664e317216ca0d687ff359a5a8be</anchor>
      <arglist>(int u_idx)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>getAdvectionVelocity</name>
      <anchorfile>class_i_b_a_m_r_1_1_convective_operator.html</anchorfile>
      <anchor>a2718b8891403591f512664eaa34a7f44</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setConvectiveDifferencingType</name>
      <anchorfile>class_i_b_a_m_r_1_1_convective_operator.html</anchorfile>
      <anchor>ae359d84ac2db09d5d51eb364cfe2819a</anchor>
      <arglist>(ConvectiveDifferencingType difference_form)</arglist>
    </member>
    <member kind="function">
      <type>ConvectiveDifferencingType</type>
      <name>getConvectiveDifferencingType</name>
      <anchorfile>class_i_b_a_m_r_1_1_convective_operator.html</anchorfile>
      <anchor>a99e0156803f41fa491f05c023d394236</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>apply</name>
      <anchorfile>class_i_b_a_m_r_1_1_convective_operator.html</anchorfile>
      <anchor>a9a031ecbac5454a8733ecac80466ff91</anchor>
      <arglist>(SAMRAI::solv::SAMRAIVectorReal&lt; NDIM, double &gt; &amp;x, SAMRAI::solv::SAMRAIVectorReal&lt; NDIM, double &gt; &amp;y)</arglist>
    </member>
    <member kind="function" static="yes">
      <type>static SAMRAI::tbox::Pointer&lt; ConvectiveOperator &gt;</type>
      <name>allocate_operator</name>
      <anchorfile>class_i_b_a_m_r_1_1_adv_diff_centered_convective_operator.html</anchorfile>
      <anchor>a4ff9bfd235e020a778084aa2a1a2f55e</anchor>
      <arglist>(const std::string &amp;object_name, SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::CellVariable&lt; NDIM, double &gt; &gt; Q_var, SAMRAI::tbox::Pointer&lt; SAMRAI::tbox::Database &gt; input_db, ConvectiveDifferencingType difference_form, const std::vector&lt; SAMRAI::solv::RobinBcCoefStrategy&lt; NDIM &gt; * &gt; &amp;bc_coefs)</arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>ConvectiveDifferencingType</type>
      <name>d_difference_form</name>
      <anchorfile>class_i_b_a_m_r_1_1_convective_operator.html</anchorfile>
      <anchor>a31a318436a5e0049d9c616b3f3756b01</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>int</type>
      <name>d_u_idx</name>
      <anchorfile>class_i_b_a_m_r_1_1_convective_operator.html</anchorfile>
      <anchor>a69bc4e145338d284a117aa106991de66</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>IBAMR::AdvDiffConvectiveOperatorManager</name>
    <filename>class_i_b_a_m_r_1_1_adv_diff_convective_operator_manager.html</filename>
    <member kind="typedef">
      <type>SAMRAI::tbox::Pointer&lt; ConvectiveOperator &gt;(*</type>
      <name>OperatorMaker</name>
      <anchorfile>class_i_b_a_m_r_1_1_adv_diff_convective_operator_manager.html</anchorfile>
      <anchor>a9ba309be369b7283e4aaa42f2292ffcb</anchor>
      <arglist>)(const std::string &amp;operator_object_name, SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::CellVariable&lt; NDIM, double &gt; &gt; Q_var, SAMRAI::tbox::Pointer&lt; SAMRAI::tbox::Database &gt; input_db, ConvectiveDifferencingType difference_form, const std::vector&lt; SAMRAI::solv::RobinBcCoefStrategy&lt; NDIM &gt; * &gt; &amp;bc_coefs)</arglist>
    </member>
    <member kind="function">
      <type>SAMRAI::tbox::Pointer&lt; ConvectiveOperator &gt;</type>
      <name>allocateOperator</name>
      <anchorfile>class_i_b_a_m_r_1_1_adv_diff_convective_operator_manager.html</anchorfile>
      <anchor>a0fc9b2eb272ec695dd988bf354a7d78d</anchor>
      <arglist>(const std::string &amp;operator_type, const std::string &amp;operator_object_name, SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::CellVariable&lt; NDIM, double &gt; &gt; Q_var, SAMRAI::tbox::Pointer&lt; SAMRAI::tbox::Database &gt; input_db, ConvectiveDifferencingType difference_form, const std::vector&lt; SAMRAI::solv::RobinBcCoefStrategy&lt; NDIM &gt; * &gt; &amp;bc_coefs) const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>registerOperatorFactoryFunction</name>
      <anchorfile>class_i_b_a_m_r_1_1_adv_diff_convective_operator_manager.html</anchorfile>
      <anchor>a615e243d40c6aa194db2171be06ca667</anchor>
      <arglist>(const std::string &amp;operator_type, OperatorMaker operator_maker)</arglist>
    </member>
    <member kind="function" static="yes">
      <type>static AdvDiffConvectiveOperatorManager *</type>
      <name>getManager</name>
      <anchorfile>class_i_b_a_m_r_1_1_adv_diff_convective_operator_manager.html</anchorfile>
      <anchor>a41c2b05621c7a9ff1fa2231788d4f1f3</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" static="yes">
      <type>static void</type>
      <name>freeManager</name>
      <anchorfile>class_i_b_a_m_r_1_1_adv_diff_convective_operator_manager.html</anchorfile>
      <anchor>a2ad6665ae1bbd7ca7788a4ed7b4a3c37</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="variable" static="yes">
      <type>static const std::string</type>
      <name>DEFAULT</name>
      <anchorfile>class_i_b_a_m_r_1_1_adv_diff_convective_operator_manager.html</anchorfile>
      <anchor>a316a0b7b51956b7eb28c9eefeaef6c5d</anchor>
      <arglist></arglist>
    </member>
    <member kind="function" protection="protected">
      <type></type>
      <name>AdvDiffConvectiveOperatorManager</name>
      <anchorfile>class_i_b_a_m_r_1_1_adv_diff_convective_operator_manager.html</anchorfile>
      <anchor>acd7732aac06cf7134c61f8adf42ac406</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" protection="protected">
      <type></type>
      <name>~AdvDiffConvectiveOperatorManager</name>
      <anchorfile>class_i_b_a_m_r_1_1_adv_diff_convective_operator_manager.html</anchorfile>
      <anchor>a4ce7638532f845bf64e3e136d25ac17c</anchor>
      <arglist>()</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>IBAMR::AdvDiffHierarchyIntegrator</name>
    <filename>class_i_b_a_m_r_1_1_adv_diff_hierarchy_integrator.html</filename>
    <base>IBTK::HierarchyIntegrator</base>
    <member kind="function">
      <type></type>
      <name>~AdvDiffHierarchyIntegrator</name>
      <anchorfile>class_i_b_a_m_r_1_1_adv_diff_hierarchy_integrator.html</anchorfile>
      <anchor>ae17d1dfd08813c1f5d6ca044e5e305a8</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setDefaultDiffusionTimeSteppingType</name>
      <anchorfile>class_i_b_a_m_r_1_1_adv_diff_hierarchy_integrator.html</anchorfile>
      <anchor>a871d656c326a20a89b93be8604b8ba68</anchor>
      <arglist>(TimeSteppingType default_diffusion_time_stepping_type)</arglist>
    </member>
    <member kind="function">
      <type>TimeSteppingType</type>
      <name>getDefaultDiffusionTimeSteppingType</name>
      <anchorfile>class_i_b_a_m_r_1_1_adv_diff_hierarchy_integrator.html</anchorfile>
      <anchor>a730a5f1ea0ef33903b8e626cd9aa3691</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setDefaultConvectiveDifferencingType</name>
      <anchorfile>class_i_b_a_m_r_1_1_adv_diff_hierarchy_integrator.html</anchorfile>
      <anchor>ac1589a985ede5413c88903cb2ace488e</anchor>
      <arglist>(ConvectiveDifferencingType default_convective_difference_form)</arglist>
    </member>
    <member kind="function">
      <type>ConvectiveDifferencingType</type>
      <name>getDefaultConvectiveDifferencingType</name>
      <anchorfile>class_i_b_a_m_r_1_1_adv_diff_hierarchy_integrator.html</anchorfile>
      <anchor>a2855c1b399b626e3ad3cbffceae96337</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>registerAdvectionVelocity</name>
      <anchorfile>class_i_b_a_m_r_1_1_adv_diff_hierarchy_integrator.html</anchorfile>
      <anchor>a6d2edaff1e9d4ffca2d9087312ca9442</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::FaceVariable&lt; NDIM, double &gt; &gt; u_var)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setAdvectionVelocityIsDivergenceFree</name>
      <anchorfile>class_i_b_a_m_r_1_1_adv_diff_hierarchy_integrator.html</anchorfile>
      <anchor>a1c3424655f8325fb00e2d148dd816465</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::FaceVariable&lt; NDIM, double &gt; &gt; u_var, bool is_div_free)</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>getAdvectionVelocityIsDivergenceFree</name>
      <anchorfile>class_i_b_a_m_r_1_1_adv_diff_hierarchy_integrator.html</anchorfile>
      <anchor>acc166e8fcbb0f6399c79368df5db35d9</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::FaceVariable&lt; NDIM, double &gt; &gt; u_var) const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setAdvectionVelocityFunction</name>
      <anchorfile>class_i_b_a_m_r_1_1_adv_diff_hierarchy_integrator.html</anchorfile>
      <anchor>a7672f90edf3a9149a3edd6979a3d9193</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::FaceVariable&lt; NDIM, double &gt; &gt; u_var, SAMRAI::tbox::Pointer&lt; IBTK::CartGridFunction &gt; u_fcn)</arglist>
    </member>
    <member kind="function">
      <type>SAMRAI::tbox::Pointer&lt; IBTK::CartGridFunction &gt;</type>
      <name>getAdvectionVelocityFunction</name>
      <anchorfile>class_i_b_a_m_r_1_1_adv_diff_hierarchy_integrator.html</anchorfile>
      <anchor>a110fc9757d06d3cffebff78264553518</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::FaceVariable&lt; NDIM, double &gt; &gt; u_var) const </arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>registerSourceTerm</name>
      <anchorfile>class_i_b_a_m_r_1_1_adv_diff_hierarchy_integrator.html</anchorfile>
      <anchor>ab43a404b44795768c5ca8b15cbdcd53c</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::CellVariable&lt; NDIM, double &gt; &gt; F_var)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setSourceTermFunction</name>
      <anchorfile>class_i_b_a_m_r_1_1_adv_diff_hierarchy_integrator.html</anchorfile>
      <anchor>aa92cd74cbabc8302c95f8006f62d339f</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::CellVariable&lt; NDIM, double &gt; &gt; F_var, SAMRAI::tbox::Pointer&lt; IBTK::CartGridFunction &gt; F_fcn)</arglist>
    </member>
    <member kind="function">
      <type>SAMRAI::tbox::Pointer&lt; IBTK::CartGridFunction &gt;</type>
      <name>getSourceTermFunction</name>
      <anchorfile>class_i_b_a_m_r_1_1_adv_diff_hierarchy_integrator.html</anchorfile>
      <anchor>a28b03b5c7a57216cec528991fa18efe6</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::CellVariable&lt; NDIM, double &gt; &gt; F_var) const </arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>registerTransportedQuantity</name>
      <anchorfile>class_i_b_a_m_r_1_1_adv_diff_hierarchy_integrator.html</anchorfile>
      <anchor>a87be97db705a8bdf3eb864c0f8cc042e</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::CellVariable&lt; NDIM, double &gt; &gt; Q_var)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setAdvectionVelocity</name>
      <anchorfile>class_i_b_a_m_r_1_1_adv_diff_hierarchy_integrator.html</anchorfile>
      <anchor>ab042f83b217d47be7e64fb77db0ca7c2</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::CellVariable&lt; NDIM, double &gt; &gt; Q_var, SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::FaceVariable&lt; NDIM, double &gt; &gt; u_var)</arglist>
    </member>
    <member kind="function">
      <type>SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::FaceVariable&lt; NDIM, double &gt; &gt;</type>
      <name>getAdvectionVelocity</name>
      <anchorfile>class_i_b_a_m_r_1_1_adv_diff_hierarchy_integrator.html</anchorfile>
      <anchor>aac56146ed4a4c370a2548ac826c98816</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::CellVariable&lt; NDIM, double &gt; &gt; Q_var) const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setSourceTerm</name>
      <anchorfile>class_i_b_a_m_r_1_1_adv_diff_hierarchy_integrator.html</anchorfile>
      <anchor>acac8f69ecdc2452794bc8d889d2183c9</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::CellVariable&lt; NDIM, double &gt; &gt; Q_var, SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::CellVariable&lt; NDIM, double &gt; &gt; F_var)</arglist>
    </member>
    <member kind="function">
      <type>SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::CellVariable&lt; NDIM, double &gt; &gt;</type>
      <name>getSourceTerm</name>
      <anchorfile>class_i_b_a_m_r_1_1_adv_diff_hierarchy_integrator.html</anchorfile>
      <anchor>ade64fcdccb00f7c7a100f5900f7c4c09</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::CellVariable&lt; NDIM, double &gt; &gt; Q_var) const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setDiffusionTimeSteppingType</name>
      <anchorfile>class_i_b_a_m_r_1_1_adv_diff_hierarchy_integrator.html</anchorfile>
      <anchor>a543de7c45c47e9bc30e193ef2185dc84</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::CellVariable&lt; NDIM, double &gt; &gt; Q_var, TimeSteppingType time_stepping_type)</arglist>
    </member>
    <member kind="function">
      <type>TimeSteppingType</type>
      <name>getDiffusionTimeSteppingType</name>
      <anchorfile>class_i_b_a_m_r_1_1_adv_diff_hierarchy_integrator.html</anchorfile>
      <anchor>a96ba2f84e3853da170ea8506742f163f</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::CellVariable&lt; NDIM, double &gt; &gt; Q_var) const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setConvectiveDifferencingType</name>
      <anchorfile>class_i_b_a_m_r_1_1_adv_diff_hierarchy_integrator.html</anchorfile>
      <anchor>ae811be5e5b8728df581705889eee5f42</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::CellVariable&lt; NDIM, double &gt; &gt; Q_var, ConvectiveDifferencingType difference_form)</arglist>
    </member>
    <member kind="function">
      <type>ConvectiveDifferencingType</type>
      <name>getConvectiveDifferencingType</name>
      <anchorfile>class_i_b_a_m_r_1_1_adv_diff_hierarchy_integrator.html</anchorfile>
      <anchor>a558e14134aed31ab4e234ab2bb3c942c</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::CellVariable&lt; NDIM, double &gt; &gt; Q_var) const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setDiffusionCoefficient</name>
      <anchorfile>class_i_b_a_m_r_1_1_adv_diff_hierarchy_integrator.html</anchorfile>
      <anchor>a4ce1aa3972125c518655141b845e37b8</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::CellVariable&lt; NDIM, double &gt; &gt; Q_var, double kappa)</arglist>
    </member>
    <member kind="function">
      <type>double</type>
      <name>getDiffusionCoefficient</name>
      <anchorfile>class_i_b_a_m_r_1_1_adv_diff_hierarchy_integrator.html</anchorfile>
      <anchor>ac0fc0534a2e4351b904f721bfd8a247e</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::CellVariable&lt; NDIM, double &gt; &gt; Q_var) const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>registerDiffusionCoefficientVariable</name>
      <anchorfile>class_i_b_a_m_r_1_1_adv_diff_hierarchy_integrator.html</anchorfile>
      <anchor>a7021c96d44e456a6a080d7a3799cdd4a</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::SideVariable&lt; NDIM, double &gt; &gt; D_var)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setDiffusionCoefficientFunction</name>
      <anchorfile>class_i_b_a_m_r_1_1_adv_diff_hierarchy_integrator.html</anchorfile>
      <anchor>a23fe83be044d4e5275163277d134b3ee</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::SideVariable&lt; NDIM, double &gt; &gt; D_var, SAMRAI::tbox::Pointer&lt; IBTK::CartGridFunction &gt; D_fcn)</arglist>
    </member>
    <member kind="function">
      <type>SAMRAI::tbox::Pointer&lt; IBTK::CartGridFunction &gt;</type>
      <name>getDiffusionCoefficientFunction</name>
      <anchorfile>class_i_b_a_m_r_1_1_adv_diff_hierarchy_integrator.html</anchorfile>
      <anchor>aac04e25fa77f6cc0e691e02f1e826741</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::SideVariable&lt; NDIM, double &gt; &gt; D_var) const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setDiffusionCoefficientVariable</name>
      <anchorfile>class_i_b_a_m_r_1_1_adv_diff_hierarchy_integrator.html</anchorfile>
      <anchor>abc6b28f3d632ce0f68a7b7b73e74d982</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::CellVariable&lt; NDIM, double &gt; &gt; Q_var, SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::SideVariable&lt; NDIM, double &gt; &gt; D_var)</arglist>
    </member>
    <member kind="function">
      <type>SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::SideVariable&lt; NDIM, double &gt; &gt;</type>
      <name>getDiffusionCoefficientVariable</name>
      <anchorfile>class_i_b_a_m_r_1_1_adv_diff_hierarchy_integrator.html</anchorfile>
      <anchor>aba3f1bffb2c024da2341733fb1091d40</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::CellVariable&lt; NDIM, double &gt; &gt; Q_var) const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setDampingCoefficient</name>
      <anchorfile>class_i_b_a_m_r_1_1_adv_diff_hierarchy_integrator.html</anchorfile>
      <anchor>ac771457764b934dfbe30d0e2c25495cf</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::CellVariable&lt; NDIM, double &gt; &gt; Q_var, double lambda)</arglist>
    </member>
    <member kind="function">
      <type>double</type>
      <name>getDampingCoefficient</name>
      <anchorfile>class_i_b_a_m_r_1_1_adv_diff_hierarchy_integrator.html</anchorfile>
      <anchor>a00257a2af962bdf97e87929a6da515f3</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::CellVariable&lt; NDIM, double &gt; &gt; Q_var) const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setInitialConditions</name>
      <anchorfile>class_i_b_a_m_r_1_1_adv_diff_hierarchy_integrator.html</anchorfile>
      <anchor>a5db49de0832de53a60248440b44378c9</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::CellVariable&lt; NDIM, double &gt; &gt; Q_var, SAMRAI::tbox::Pointer&lt; IBTK::CartGridFunction &gt; Q_init)</arglist>
    </member>
    <member kind="function">
      <type>SAMRAI::tbox::Pointer&lt; IBTK::CartGridFunction &gt;</type>
      <name>getInitialConditions</name>
      <anchorfile>class_i_b_a_m_r_1_1_adv_diff_hierarchy_integrator.html</anchorfile>
      <anchor>ae93c7f2623ea40a6d6fb230004f2d6a9</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::CellVariable&lt; NDIM, double &gt; &gt; Q_var) const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setPhysicalBcCoef</name>
      <anchorfile>class_i_b_a_m_r_1_1_adv_diff_hierarchy_integrator.html</anchorfile>
      <anchor>ab6e6bfee80bfc1fa15cfa0fb1399889b</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::CellVariable&lt; NDIM, double &gt; &gt; Q_var, SAMRAI::solv::RobinBcCoefStrategy&lt; NDIM &gt; *Q_bc_coef)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setPhysicalBcCoefs</name>
      <anchorfile>class_i_b_a_m_r_1_1_adv_diff_hierarchy_integrator.html</anchorfile>
      <anchor>a0735e067b2bfebd348b15074b0df259a</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::CellVariable&lt; NDIM, double &gt; &gt; Q_var, const std::vector&lt; SAMRAI::solv::RobinBcCoefStrategy&lt; NDIM &gt; * &gt; &amp;Q_bc_coef)</arglist>
    </member>
    <member kind="function">
      <type>std::vector&lt; SAMRAI::solv::RobinBcCoefStrategy&lt; NDIM &gt; * &gt;</type>
      <name>getPhysicalBcCoefs</name>
      <anchorfile>class_i_b_a_m_r_1_1_adv_diff_hierarchy_integrator.html</anchorfile>
      <anchor>a84fbf71e866cb32fee4a6f89a9470676</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::CellVariable&lt; NDIM, double &gt; &gt; Q_var) const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setHelmholtzSolver</name>
      <anchorfile>class_i_b_a_m_r_1_1_adv_diff_hierarchy_integrator.html</anchorfile>
      <anchor>aac933584f985ec798d8700a56119de67</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::CellVariable&lt; NDIM, double &gt; &gt; Q_var, SAMRAI::tbox::Pointer&lt; IBTK::PoissonSolver &gt; helmholtz_solver)</arglist>
    </member>
    <member kind="function">
      <type>SAMRAI::tbox::Pointer&lt; IBTK::PoissonSolver &gt;</type>
      <name>getHelmholtzSolver</name>
      <anchorfile>class_i_b_a_m_r_1_1_adv_diff_hierarchy_integrator.html</anchorfile>
      <anchor>a6280bc438e8638cd169e8bf01bb7cf51</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::CellVariable&lt; NDIM, double &gt; &gt; Q_var)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setHelmholtzRHSOperator</name>
      <anchorfile>class_i_b_a_m_r_1_1_adv_diff_hierarchy_integrator.html</anchorfile>
      <anchor>ac9c18e38257f11be4029a775b7e2d00d</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::CellVariable&lt; NDIM, double &gt; &gt; Q_var, SAMRAI::tbox::Pointer&lt; IBTK::LaplaceOperator &gt; helmholtz_rhs_operator)</arglist>
    </member>
    <member kind="function">
      <type>SAMRAI::tbox::Pointer&lt; IBTK::LaplaceOperator &gt;</type>
      <name>getHelmholtzRHSOperator</name>
      <anchorfile>class_i_b_a_m_r_1_1_adv_diff_hierarchy_integrator.html</anchorfile>
      <anchor>aeb8a86926fc7680e9ae666b2fd2e0653</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::CellVariable&lt; NDIM, double &gt; &gt; Q_var)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>initializeHierarchyIntegrator</name>
      <anchorfile>class_i_b_a_m_r_1_1_adv_diff_hierarchy_integrator.html</anchorfile>
      <anchor>af7c27b14476594357049ccbd28b2fd7a</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::hier::PatchHierarchy&lt; NDIM &gt; &gt; hierarchy, SAMRAI::tbox::Pointer&lt; SAMRAI::mesh::GriddingAlgorithm&lt; NDIM &gt; &gt; gridding_alg)</arglist>
    </member>
    <member kind="function" protection="protected">
      <type></type>
      <name>AdvDiffHierarchyIntegrator</name>
      <anchorfile>class_i_b_a_m_r_1_1_adv_diff_hierarchy_integrator.html</anchorfile>
      <anchor>a2683c6ab9c73d22813707f59dbbff8f7</anchor>
      <arglist>(const std::string &amp;object_name, SAMRAI::tbox::Pointer&lt; SAMRAI::tbox::Database &gt; input_db, bool register_for_restart)</arglist>
    </member>
    <member kind="function" protection="protected">
      <type>double</type>
      <name>getMaximumTimeStepSizeSpecialized</name>
      <anchorfile>class_i_b_a_m_r_1_1_adv_diff_hierarchy_integrator.html</anchorfile>
      <anchor>a4f64c2b6c8d2a691448bcb148353cd0c</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" protection="protected">
      <type>void</type>
      <name>resetHierarchyConfigurationSpecialized</name>
      <anchorfile>class_i_b_a_m_r_1_1_adv_diff_hierarchy_integrator.html</anchorfile>
      <anchor>a8657e002e672f619052f495b32e6beef</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::hier::BasePatchHierarchy&lt; NDIM &gt; &gt; hierarchy, int coarsest_level, int finest_level)</arglist>
    </member>
    <member kind="function" protection="protected">
      <type>void</type>
      <name>putToDatabaseSpecialized</name>
      <anchorfile>class_i_b_a_m_r_1_1_adv_diff_hierarchy_integrator.html</anchorfile>
      <anchor>af9eb41fe4b7d2c3217dbaf85376bb870</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::tbox::Database &gt; db)</arglist>
    </member>
    <member kind="function" protection="protected">
      <type>void</type>
      <name>registerVariables</name>
      <anchorfile>class_i_b_a_m_r_1_1_adv_diff_hierarchy_integrator.html</anchorfile>
      <anchor>af2e495595ba57e05ecb999a35fbf23d8</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>double</type>
      <name>d_cfl_max</name>
      <anchorfile>class_i_b_a_m_r_1_1_adv_diff_hierarchy_integrator.html</anchorfile>
      <anchor>ad6bcccd9f37f64bad79695e03aa61ca2</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>TimeSteppingType</type>
      <name>d_default_diffusion_time_stepping_type</name>
      <anchorfile>class_i_b_a_m_r_1_1_adv_diff_hierarchy_integrator.html</anchorfile>
      <anchor>a4cc59aef746d6847318aa55a3c7cc640</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>ConvectiveDifferencingType</type>
      <name>d_default_convective_difference_form</name>
      <anchorfile>class_i_b_a_m_r_1_1_adv_diff_hierarchy_integrator.html</anchorfile>
      <anchor>a668b64debc95481be87ee155b1827cfa</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>std::vector&lt; SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::FaceVariable&lt; NDIM, double &gt; &gt; &gt;</type>
      <name>d_u_var</name>
      <anchorfile>class_i_b_a_m_r_1_1_adv_diff_hierarchy_integrator.html</anchorfile>
      <anchor>a4e5fb53710a42721e40fb9edb841add1</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>std::vector&lt; SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::CellVariable&lt; NDIM, double &gt; &gt; &gt;</type>
      <name>d_F_var</name>
      <anchorfile>class_i_b_a_m_r_1_1_adv_diff_hierarchy_integrator.html</anchorfile>
      <anchor>a5a6786f843dc0089052e72ab7dab7b15</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>std::vector&lt; SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::SideVariable&lt; NDIM, double &gt; &gt; &gt;</type>
      <name>d_diffusion_coef_var</name>
      <anchorfile>class_i_b_a_m_r_1_1_adv_diff_hierarchy_integrator.html</anchorfile>
      <anchor>afe589dedd1821c0106ea74ca27e84ffc</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>std::vector&lt; SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::CellVariable&lt; NDIM, double &gt; &gt; &gt;</type>
      <name>d_Q_var</name>
      <anchorfile>class_i_b_a_m_r_1_1_adv_diff_hierarchy_integrator.html</anchorfile>
      <anchor>af1e8945099bc4fa88951a6ece189377b</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>IBAMR::AdvDiffPhysicalBoundaryUtilities</name>
    <filename>class_i_b_a_m_r_1_1_adv_diff_physical_boundary_utilities.html</filename>
    <member kind="function" static="yes">
      <type>static void</type>
      <name>setPhysicalBoundaryConditions</name>
      <anchorfile>class_i_b_a_m_r_1_1_adv_diff_physical_boundary_utilities.html</anchorfile>
      <anchor>a146753b5d9232893f4d4d567550ff046</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::CellData&lt; NDIM, double &gt; &gt; Q_data, SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::FaceData&lt; NDIM, double &gt; &gt; u_ADV_data, SAMRAI::tbox::Pointer&lt; SAMRAI::hier::Patch&lt; NDIM &gt; &gt; patch, const std::vector&lt; SAMRAI::solv::RobinBcCoefStrategy&lt; NDIM &gt; * &gt; &amp;bc_coefs, const double fill_time, const bool inflow_boundaries_only, const bool homogeneous_bc)</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>IBAMR::AdvDiffPPMConvectiveOperator</name>
    <filename>class_i_b_a_m_r_1_1_adv_diff_p_p_m_convective_operator.html</filename>
    <base>IBAMR::ConvectiveOperator</base>
    <member kind="function">
      <type></type>
      <name>AdvDiffPPMConvectiveOperator</name>
      <anchorfile>class_i_b_a_m_r_1_1_adv_diff_p_p_m_convective_operator.html</anchorfile>
      <anchor>ab5d0a3f6c30894c0c439e9ba0788d7d5</anchor>
      <arglist>(const std::string &amp;object_name, SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::CellVariable&lt; NDIM, double &gt; &gt; Q_var, SAMRAI::tbox::Pointer&lt; SAMRAI::tbox::Database &gt; input_db, ConvectiveDifferencingType difference_form, const std::vector&lt; SAMRAI::solv::RobinBcCoefStrategy&lt; NDIM &gt; * &gt; &amp;bc_coefs)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~AdvDiffPPMConvectiveOperator</name>
      <anchorfile>class_i_b_a_m_r_1_1_adv_diff_p_p_m_convective_operator.html</anchorfile>
      <anchor>ae99fb2aa7a2a15309755bafb1f44fb13</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>applyConvectiveOperator</name>
      <anchorfile>class_i_b_a_m_r_1_1_adv_diff_p_p_m_convective_operator.html</anchorfile>
      <anchor>a43284842f57fe1208a3a9cee47344104</anchor>
      <arglist>(int Q_idx, int N_idx)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>initializeOperatorState</name>
      <anchorfile>class_i_b_a_m_r_1_1_adv_diff_p_p_m_convective_operator.html</anchorfile>
      <anchor>a68c33d733a6faebcd3605588d8c1c62e</anchor>
      <arglist>(const SAMRAI::solv::SAMRAIVectorReal&lt; NDIM, double &gt; &amp;in, const SAMRAI::solv::SAMRAIVectorReal&lt; NDIM, double &gt; &amp;out)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>deallocateOperatorState</name>
      <anchorfile>class_i_b_a_m_r_1_1_adv_diff_p_p_m_convective_operator.html</anchorfile>
      <anchor>ad619cf6234e950f02826f8f038442149</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" static="yes">
      <type>static SAMRAI::tbox::Pointer&lt; ConvectiveOperator &gt;</type>
      <name>allocate_operator</name>
      <anchorfile>class_i_b_a_m_r_1_1_adv_diff_p_p_m_convective_operator.html</anchorfile>
      <anchor>a86cf21afb11e3af6148fea58da05bdae</anchor>
      <arglist>(const std::string &amp;object_name, SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::CellVariable&lt; NDIM, double &gt; &gt; Q_var, SAMRAI::tbox::Pointer&lt; SAMRAI::tbox::Database &gt; input_db, ConvectiveDifferencingType difference_form, const std::vector&lt; SAMRAI::solv::RobinBcCoefStrategy&lt; NDIM &gt; * &gt; &amp;bc_coefs)</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>IBAMR::AdvDiffPredictorCorrectorHierarchyIntegrator</name>
    <filename>class_i_b_a_m_r_1_1_adv_diff_predictor_corrector_hierarchy_integrator.html</filename>
    <base>IBAMR::AdvDiffHierarchyIntegrator</base>
    <member kind="function">
      <type></type>
      <name>AdvDiffPredictorCorrectorHierarchyIntegrator</name>
      <anchorfile>class_i_b_a_m_r_1_1_adv_diff_predictor_corrector_hierarchy_integrator.html</anchorfile>
      <anchor>ae2cfcc70ad306e8044976afd3f3edbdf</anchor>
      <arglist>(const std::string &amp;object_name, SAMRAI::tbox::Pointer&lt; SAMRAI::tbox::Database &gt; input_db, SAMRAI::tbox::Pointer&lt; AdvectorExplicitPredictorPatchOps &gt; explicit_predictor, bool register_for_restart=true)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~AdvDiffPredictorCorrectorHierarchyIntegrator</name>
      <anchorfile>class_i_b_a_m_r_1_1_adv_diff_predictor_corrector_hierarchy_integrator.html</anchorfile>
      <anchor>a428cab98fadb14f96cc3cf9101b034a7</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>SAMRAI::tbox::Pointer&lt; SAMRAI::algs::HyperbolicLevelIntegrator&lt; NDIM &gt; &gt;</type>
      <name>getHyperbolicLevelIntegrator</name>
      <anchorfile>class_i_b_a_m_r_1_1_adv_diff_predictor_corrector_hierarchy_integrator.html</anchorfile>
      <anchor>a9cda4d2b62e97bf4be0b8b277b98e381</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>SAMRAI::tbox::Pointer&lt; AdvDiffPredictorCorrectorHyperbolicPatchOps &gt;</type>
      <name>getHyperbolicPatchStrategy</name>
      <anchorfile>class_i_b_a_m_r_1_1_adv_diff_predictor_corrector_hierarchy_integrator.html</anchorfile>
      <anchor>a83275a796b8475bf73999cf6053b749e</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>initializeHierarchyIntegrator</name>
      <anchorfile>class_i_b_a_m_r_1_1_adv_diff_predictor_corrector_hierarchy_integrator.html</anchorfile>
      <anchor>a9678c741e5b428e44c46e5a550e6cb7d</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::hier::PatchHierarchy&lt; NDIM &gt; &gt; hierarchy, SAMRAI::tbox::Pointer&lt; SAMRAI::mesh::GriddingAlgorithm&lt; NDIM &gt; &gt; gridding_alg)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>preprocessIntegrateHierarchy</name>
      <anchorfile>class_i_b_a_m_r_1_1_adv_diff_predictor_corrector_hierarchy_integrator.html</anchorfile>
      <anchor>abdc1ca6ae6c6305dfed6509ff3b2f402</anchor>
      <arglist>(double current_time, double new_time, int num_cycles=1)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>integrateHierarchy</name>
      <anchorfile>class_i_b_a_m_r_1_1_adv_diff_predictor_corrector_hierarchy_integrator.html</anchorfile>
      <anchor>a34939f77e8d8cf6ee4d33bdad653e49d</anchor>
      <arglist>(double current_time, double new_time, int cycle_num=0)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>postprocessIntegrateHierarchy</name>
      <anchorfile>class_i_b_a_m_r_1_1_adv_diff_predictor_corrector_hierarchy_integrator.html</anchorfile>
      <anchor>a7c98e1acf701833179423ca0844f4329</anchor>
      <arglist>(double current_time, double new_time, bool skip_synchronize_new_state_data, int num_cycles=1)</arglist>
    </member>
    <member kind="function" protection="protected">
      <type>double</type>
      <name>getMaximumTimeStepSizeSpecialized</name>
      <anchorfile>class_i_b_a_m_r_1_1_adv_diff_predictor_corrector_hierarchy_integrator.html</anchorfile>
      <anchor>a0434dddc72b8b25eee2b03e5af2a0bcf</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" protection="protected">
      <type>void</type>
      <name>resetTimeDependentHierarchyDataSpecialized</name>
      <anchorfile>class_i_b_a_m_r_1_1_adv_diff_predictor_corrector_hierarchy_integrator.html</anchorfile>
      <anchor>a50b4098aadec5f6adbd97bfb3e9b6285</anchor>
      <arglist>(double new_time)</arglist>
    </member>
    <member kind="function" protection="protected">
      <type>void</type>
      <name>resetIntegratorToPreadvanceStateSpecialized</name>
      <anchorfile>class_i_b_a_m_r_1_1_adv_diff_predictor_corrector_hierarchy_integrator.html</anchorfile>
      <anchor>a6dadf0898c6b9a7b610424345edc25d7</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" protection="protected">
      <type>void</type>
      <name>initializeLevelDataSpecialized</name>
      <anchorfile>class_i_b_a_m_r_1_1_adv_diff_predictor_corrector_hierarchy_integrator.html</anchorfile>
      <anchor>aee596677415de525cf318c9bd5371ee4</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::hier::BasePatchHierarchy&lt; NDIM &gt; &gt; hierarchy, int level_number, double init_data_time, bool can_be_refined, bool initial_time, SAMRAI::tbox::Pointer&lt; SAMRAI::hier::BasePatchLevel&lt; NDIM &gt; &gt; old_level, bool allocate_data)</arglist>
    </member>
    <member kind="function" protection="protected">
      <type>void</type>
      <name>resetHierarchyConfigurationSpecialized</name>
      <anchorfile>class_i_b_a_m_r_1_1_adv_diff_predictor_corrector_hierarchy_integrator.html</anchorfile>
      <anchor>ade3aef98da3f3a4528e0a829a8f914c5</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::hier::BasePatchHierarchy&lt; NDIM &gt; &gt; hierarchy, int coarsest_level, int finest_level)</arglist>
    </member>
    <member kind="function" protection="protected">
      <type>void</type>
      <name>applyGradientDetectorSpecialized</name>
      <anchorfile>class_i_b_a_m_r_1_1_adv_diff_predictor_corrector_hierarchy_integrator.html</anchorfile>
      <anchor>ad0ea832cdf14fe3f1344dd1dd3e8adfa</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::hier::BasePatchHierarchy&lt; NDIM &gt; &gt; hierarchy, int level_number, double error_data_time, int tag_index, bool initial_time, bool uses_richardson_extrapolation_too)</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>IBAMR::AdvDiffPredictorCorrectorHyperbolicPatchOps</name>
    <filename>class_i_b_a_m_r_1_1_adv_diff_predictor_corrector_hyperbolic_patch_ops.html</filename>
    <base>IBAMR::AdvectorPredictorCorrectorHyperbolicPatchOps</base>
    <member kind="function">
      <type></type>
      <name>AdvDiffPredictorCorrectorHyperbolicPatchOps</name>
      <anchorfile>class_i_b_a_m_r_1_1_adv_diff_predictor_corrector_hyperbolic_patch_ops.html</anchorfile>
      <anchor>a4b787567f9675d1ae11bcf8e7800db29</anchor>
      <arglist>(const std::string &amp;object_name, SAMRAI::tbox::Pointer&lt; SAMRAI::tbox::Database &gt; input_db, SAMRAI::tbox::Pointer&lt; AdvectorExplicitPredictorPatchOps &gt; explicit_predictor, SAMRAI::tbox::Pointer&lt; SAMRAI::geom::CartesianGridGeometry&lt; NDIM &gt; &gt; grid_geom, bool register_for_restart=true)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~AdvDiffPredictorCorrectorHyperbolicPatchOps</name>
      <anchorfile>class_i_b_a_m_r_1_1_adv_diff_predictor_corrector_hyperbolic_patch_ops.html</anchorfile>
      <anchor>a1066574f0173193659e0b45b955cfc5c</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>conservativeDifferenceOnPatch</name>
      <anchorfile>class_i_b_a_m_r_1_1_adv_diff_predictor_corrector_hyperbolic_patch_ops.html</anchorfile>
      <anchor>ae7f9904b9ae9e491c573bda994e37e33</anchor>
      <arglist>(SAMRAI::hier::Patch&lt; NDIM &gt; &amp;patch, double time, double dt, bool at_synchronization)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>preprocessAdvanceLevelState</name>
      <anchorfile>class_i_b_a_m_r_1_1_adv_diff_predictor_corrector_hyperbolic_patch_ops.html</anchorfile>
      <anchor>a359856474616cc8082e2d004aafd26d0</anchor>
      <arglist>(const SAMRAI::tbox::Pointer&lt; SAMRAI::hier::PatchLevel&lt; NDIM &gt; &gt; &amp;level, double current_time, double dt, bool first_step, bool last_step, bool regrid_advance)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>postprocessAdvanceLevelState</name>
      <anchorfile>class_i_b_a_m_r_1_1_adv_diff_predictor_corrector_hyperbolic_patch_ops.html</anchorfile>
      <anchor>aea93bd5790baa9501cb5655d717ea187</anchor>
      <arglist>(const SAMRAI::tbox::Pointer&lt; SAMRAI::hier::PatchLevel&lt; NDIM &gt; &gt; &amp;level, double current_time, double dt, bool first_step, bool last_step, bool regrid_advance)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>AdvectorPredictorCorrectorHyperbolicPatchOps</name>
      <anchorfile>class_i_b_a_m_r_1_1_advector_predictor_corrector_hyperbolic_patch_ops.html</anchorfile>
      <anchor>a0111355e8ea9f9840e33f0f38a1fa412</anchor>
      <arglist>(const std::string &amp;object_name, SAMRAI::tbox::Pointer&lt; SAMRAI::tbox::Database &gt; input_db, SAMRAI::tbox::Pointer&lt; AdvectorExplicitPredictorPatchOps &gt; explicit_predictor, SAMRAI::tbox::Pointer&lt; SAMRAI::geom::CartesianGridGeometry&lt; NDIM &gt; &gt; grid_geom, bool register_for_restart=true)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual</type>
      <name>~AdvectorPredictorCorrectorHyperbolicPatchOps</name>
      <anchorfile>class_i_b_a_m_r_1_1_advector_predictor_corrector_hyperbolic_patch_ops.html</anchorfile>
      <anchor>a307fb26bef2e16684d86b0c5588455d3</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>const std::string &amp;</type>
      <name>getName</name>
      <anchorfile>class_i_b_a_m_r_1_1_advector_predictor_corrector_hyperbolic_patch_ops.html</anchorfile>
      <anchor>a8ce072a46058a1d5dfd23e66f04a85be</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>registerVisItDataWriter</name>
      <anchorfile>class_i_b_a_m_r_1_1_advector_predictor_corrector_hyperbolic_patch_ops.html</anchorfile>
      <anchor>a04225588a61746ea5c68486f66e7e5c1</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::appu::VisItDataWriter&lt; NDIM &gt; &gt; visit_writer)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>registerAdvectionVelocity</name>
      <anchorfile>class_i_b_a_m_r_1_1_advector_predictor_corrector_hyperbolic_patch_ops.html</anchorfile>
      <anchor>a257bd36364f6ec02e6e555162da6ee15</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::FaceVariable&lt; NDIM, double &gt; &gt; u_var)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setAdvectionVelocityIsDivergenceFree</name>
      <anchorfile>class_i_b_a_m_r_1_1_advector_predictor_corrector_hyperbolic_patch_ops.html</anchorfile>
      <anchor>ac613a31f56ed6af36944c42252446ee9</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::FaceVariable&lt; NDIM, double &gt; &gt; u_var, bool is_div_free)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setAdvectionVelocityFunction</name>
      <anchorfile>class_i_b_a_m_r_1_1_advector_predictor_corrector_hyperbolic_patch_ops.html</anchorfile>
      <anchor>a06f2a1c7ba81b3b9762cecec7f8e19a5</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::FaceVariable&lt; NDIM, double &gt; &gt; u_var, SAMRAI::tbox::Pointer&lt; IBTK::CartGridFunction &gt; u_fcn)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>registerSourceTerm</name>
      <anchorfile>class_i_b_a_m_r_1_1_advector_predictor_corrector_hyperbolic_patch_ops.html</anchorfile>
      <anchor>ae488d4685d767b1144c56012839f3068</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::CellVariable&lt; NDIM, double &gt; &gt; F_var)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setSourceTermFunction</name>
      <anchorfile>class_i_b_a_m_r_1_1_advector_predictor_corrector_hyperbolic_patch_ops.html</anchorfile>
      <anchor>af23651ffabb8d7b69e936dcfb22d8c13</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::CellVariable&lt; NDIM, double &gt; &gt; F_var, SAMRAI::tbox::Pointer&lt; IBTK::CartGridFunction &gt; F_fcn)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>registerTransportedQuantity</name>
      <anchorfile>class_i_b_a_m_r_1_1_advector_predictor_corrector_hyperbolic_patch_ops.html</anchorfile>
      <anchor>ad813acd9b9915cbf3e6bb10e4db5e182</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::CellVariable&lt; NDIM, double &gt; &gt; Q_var)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setAdvectionVelocity</name>
      <anchorfile>class_i_b_a_m_r_1_1_advector_predictor_corrector_hyperbolic_patch_ops.html</anchorfile>
      <anchor>a070aaf26f459edf59eac9f345a0fac4e</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::CellVariable&lt; NDIM, double &gt; &gt; Q_var, SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::FaceVariable&lt; NDIM, double &gt; &gt; u_var)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setSourceTerm</name>
      <anchorfile>class_i_b_a_m_r_1_1_advector_predictor_corrector_hyperbolic_patch_ops.html</anchorfile>
      <anchor>a34abe9bad0a7ce45a12be66b30eb153d</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::CellVariable&lt; NDIM, double &gt; &gt; Q_var, SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::CellVariable&lt; NDIM, double &gt; &gt; F_var)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setConvectiveDifferencingType</name>
      <anchorfile>class_i_b_a_m_r_1_1_advector_predictor_corrector_hyperbolic_patch_ops.html</anchorfile>
      <anchor>ae393d037104f629c428794074f93ffdb</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::CellVariable&lt; NDIM, double &gt; &gt; Q_var, ConvectiveDifferencingType difference_form)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setInitialConditions</name>
      <anchorfile>class_i_b_a_m_r_1_1_advector_predictor_corrector_hyperbolic_patch_ops.html</anchorfile>
      <anchor>ad5b40b81b49d6dc791eebc812973d6a6</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::CellVariable&lt; NDIM, double &gt; &gt; Q_var, SAMRAI::tbox::Pointer&lt; IBTK::CartGridFunction &gt; Q_init)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setPhysicalBcCoefs</name>
      <anchorfile>class_i_b_a_m_r_1_1_advector_predictor_corrector_hyperbolic_patch_ops.html</anchorfile>
      <anchor>a93956786dec639c102cd20c9d5abf830</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::CellVariable&lt; NDIM, double &gt; &gt; Q_var, SAMRAI::solv::RobinBcCoefStrategy&lt; NDIM &gt; *Q_bc_coef)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setPhysicalBcCoefs</name>
      <anchorfile>class_i_b_a_m_r_1_1_advector_predictor_corrector_hyperbolic_patch_ops.html</anchorfile>
      <anchor>a37713cb046244d6c6d070505acb36c46</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::CellVariable&lt; NDIM, double &gt; &gt; Q_var, std::vector&lt; SAMRAI::solv::RobinBcCoefStrategy&lt; NDIM &gt; * &gt; Q_bc_coef)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>registerModelVariables</name>
      <anchorfile>class_i_b_a_m_r_1_1_advector_predictor_corrector_hyperbolic_patch_ops.html</anchorfile>
      <anchor>a7c01907579be6e5e70095c33ae004f68</anchor>
      <arglist>(SAMRAI::algs::HyperbolicLevelIntegrator&lt; NDIM &gt; *integrator)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>initializeDataOnPatch</name>
      <anchorfile>class_i_b_a_m_r_1_1_advector_predictor_corrector_hyperbolic_patch_ops.html</anchorfile>
      <anchor>a0eff48bc3b59f9f40e4c6a6f25a9c654</anchor>
      <arglist>(SAMRAI::hier::Patch&lt; NDIM &gt; &amp;patch, double data_time, bool initial_time)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual double</type>
      <name>computeStableDtOnPatch</name>
      <anchorfile>class_i_b_a_m_r_1_1_advector_predictor_corrector_hyperbolic_patch_ops.html</anchorfile>
      <anchor>a0137a0c28960d44d4870a837e56ff21d</anchor>
      <arglist>(SAMRAI::hier::Patch&lt; NDIM &gt; &amp;patch, bool initial_time, double dt_time)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>computeFluxesOnPatch</name>
      <anchorfile>class_i_b_a_m_r_1_1_advector_predictor_corrector_hyperbolic_patch_ops.html</anchorfile>
      <anchor>a7a748a5baf92dcd6ffa113e91f7569a2</anchor>
      <arglist>(SAMRAI::hier::Patch&lt; NDIM &gt; &amp;patch, double time, double dt)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>tagGradientDetectorCells</name>
      <anchorfile>class_i_b_a_m_r_1_1_advector_predictor_corrector_hyperbolic_patch_ops.html</anchorfile>
      <anchor>af9a51e25115edefcee197ac9c939f8e1</anchor>
      <arglist>(SAMRAI::hier::Patch&lt; NDIM &gt; &amp;patch, double regrid_time, bool initial_error, int tag_indexindx, bool uses_richardson_extrapolation_too)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>setPhysicalBoundaryConditions</name>
      <anchorfile>class_i_b_a_m_r_1_1_advector_predictor_corrector_hyperbolic_patch_ops.html</anchorfile>
      <anchor>af7e84fe18692baca2407a9f2a8eb84ca</anchor>
      <arglist>(SAMRAI::hier::Patch&lt; NDIM &gt; &amp;patch, double fill_time, const SAMRAI::hier::IntVector&lt; NDIM &gt; &amp;ghost_width_to_fill)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>putToDatabase</name>
      <anchorfile>class_i_b_a_m_r_1_1_advector_predictor_corrector_hyperbolic_patch_ops.html</anchorfile>
      <anchor>aeee98d6154d50bd26e5ddfcad2ed2a85</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::tbox::Database &gt; db)</arglist>
    </member>
    <member kind="function" protection="protected">
      <type>SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::FaceData&lt; NDIM, double &gt; &gt;</type>
      <name>getFluxIntegralData</name>
      <anchorfile>class_i_b_a_m_r_1_1_advector_predictor_corrector_hyperbolic_patch_ops.html</anchorfile>
      <anchor>a4fad3b0a99ce62231b3fb17d0261675b</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::CellVariable&lt; NDIM, double &gt; &gt; Q_var, SAMRAI::hier::Patch&lt; NDIM &gt; &amp;patch, SAMRAI::tbox::Pointer&lt; SAMRAI::hier::VariableContext &gt; context)</arglist>
    </member>
    <member kind="function" protection="protected">
      <type>SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::FaceData&lt; NDIM, double &gt; &gt;</type>
      <name>getQIntegralData</name>
      <anchorfile>class_i_b_a_m_r_1_1_advector_predictor_corrector_hyperbolic_patch_ops.html</anchorfile>
      <anchor>ad1110df212d46be26cc977b7523067ac</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::CellVariable&lt; NDIM, double &gt; &gt; Q_var, SAMRAI::hier::Patch&lt; NDIM &gt; &amp;patch, SAMRAI::tbox::Pointer&lt; SAMRAI::hier::VariableContext &gt; context)</arglist>
    </member>
    <member kind="function" protection="protected">
      <type>SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::FaceData&lt; NDIM, double &gt; &gt;</type>
      <name>getUIntegralData</name>
      <anchorfile>class_i_b_a_m_r_1_1_advector_predictor_corrector_hyperbolic_patch_ops.html</anchorfile>
      <anchor>a7785af6147668976ca5cbcef283894b2</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::CellVariable&lt; NDIM, double &gt; &gt; Q_var, SAMRAI::hier::Patch&lt; NDIM &gt; &amp;patch, SAMRAI::tbox::Pointer&lt; SAMRAI::hier::VariableContext &gt; context)</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>IBAMR::AdvDiffSemiImplicitHierarchyIntegrator</name>
    <filename>class_i_b_a_m_r_1_1_adv_diff_semi_implicit_hierarchy_integrator.html</filename>
    <base>IBAMR::AdvDiffHierarchyIntegrator</base>
    <member kind="function">
      <type></type>
      <name>AdvDiffSemiImplicitHierarchyIntegrator</name>
      <anchorfile>class_i_b_a_m_r_1_1_adv_diff_semi_implicit_hierarchy_integrator.html</anchorfile>
      <anchor>a9146d7ed01f1e1f4b6aca3237593419e</anchor>
      <arglist>(const std::string &amp;object_name, SAMRAI::tbox::Pointer&lt; SAMRAI::tbox::Database &gt; input_db, bool register_for_restart=true)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~AdvDiffSemiImplicitHierarchyIntegrator</name>
      <anchorfile>class_i_b_a_m_r_1_1_adv_diff_semi_implicit_hierarchy_integrator.html</anchorfile>
      <anchor>a63e3ac6768968880a8a8da76c97f46bc</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setDefaultConvectiveTimeSteppingType</name>
      <anchorfile>class_i_b_a_m_r_1_1_adv_diff_semi_implicit_hierarchy_integrator.html</anchorfile>
      <anchor>aff02cdfafaab5ddbfe4c7c4679ff73c9</anchor>
      <arglist>(TimeSteppingType default_convective_time_stepping_type)</arglist>
    </member>
    <member kind="function">
      <type>TimeSteppingType</type>
      <name>getDefaultConvectiveTimeSteppingType</name>
      <anchorfile>class_i_b_a_m_r_1_1_adv_diff_semi_implicit_hierarchy_integrator.html</anchorfile>
      <anchor>a691103ea97c2d678f1a536a514c9234f</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setDefaultInitialConvectiveTimeSteppingType</name>
      <anchorfile>class_i_b_a_m_r_1_1_adv_diff_semi_implicit_hierarchy_integrator.html</anchorfile>
      <anchor>a2ff3bf4839f635ff23c5360f5c8379f5</anchor>
      <arglist>(TimeSteppingType default_init_convective_time_stepping_type)</arglist>
    </member>
    <member kind="function">
      <type>TimeSteppingType</type>
      <name>getDefaultInitialConvectiveTimeSteppingType</name>
      <anchorfile>class_i_b_a_m_r_1_1_adv_diff_semi_implicit_hierarchy_integrator.html</anchorfile>
      <anchor>a44a7fc7677fb181f68e6e890b6e2bd7c</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setDefaultConvectiveOperatorType</name>
      <anchorfile>class_i_b_a_m_r_1_1_adv_diff_semi_implicit_hierarchy_integrator.html</anchorfile>
      <anchor>ab8baba7df2ea28a847b7c5f13e23cceb</anchor>
      <arglist>(const std::string &amp;op_type)</arglist>
    </member>
    <member kind="function">
      <type>const std::string &amp;</type>
      <name>getDefaultConvectiveOperatorType</name>
      <anchorfile>class_i_b_a_m_r_1_1_adv_diff_semi_implicit_hierarchy_integrator.html</anchorfile>
      <anchor>aa7f8c537a0f4ed0cf2a16a4b496acfac</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setDefaultConvectiveOperatorInputDatabase</name>
      <anchorfile>class_i_b_a_m_r_1_1_adv_diff_semi_implicit_hierarchy_integrator.html</anchorfile>
      <anchor>aa27c6432e95f9cf0ed161bdcd056b199</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::tbox::Database &gt; input_db)</arglist>
    </member>
    <member kind="function">
      <type>SAMRAI::tbox::Pointer&lt; SAMRAI::tbox::Database &gt;</type>
      <name>getDefaultConvectiveOperatorInputDatabase</name>
      <anchorfile>class_i_b_a_m_r_1_1_adv_diff_semi_implicit_hierarchy_integrator.html</anchorfile>
      <anchor>a36cfc2a9383de72de3ea01207981d28a</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>registerTransportedQuantity</name>
      <anchorfile>class_i_b_a_m_r_1_1_adv_diff_semi_implicit_hierarchy_integrator.html</anchorfile>
      <anchor>a7eccaca64f536a15d8e813638eaa83e5</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::CellVariable&lt; NDIM, double &gt; &gt; Q_var)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setConvectiveTimeSteppingType</name>
      <anchorfile>class_i_b_a_m_r_1_1_adv_diff_semi_implicit_hierarchy_integrator.html</anchorfile>
      <anchor>a8077651ec664ae2c11371e8bb427cb54</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::CellVariable&lt; NDIM, double &gt; &gt; Q_var, TimeSteppingType convective_time_stepping_type)</arglist>
    </member>
    <member kind="function">
      <type>TimeSteppingType</type>
      <name>getConvectiveTimeSteppingType</name>
      <anchorfile>class_i_b_a_m_r_1_1_adv_diff_semi_implicit_hierarchy_integrator.html</anchorfile>
      <anchor>a6e6d1e595a6153f3f87c28f091139658</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::CellVariable&lt; NDIM, double &gt; &gt; Q_var) const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setInitialConvectiveTimeSteppingType</name>
      <anchorfile>class_i_b_a_m_r_1_1_adv_diff_semi_implicit_hierarchy_integrator.html</anchorfile>
      <anchor>a4e0f53be637331012b36cdfc109d9359</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::CellVariable&lt; NDIM, double &gt; &gt; Q_var, TimeSteppingType init_convective_time_stepping_type)</arglist>
    </member>
    <member kind="function">
      <type>TimeSteppingType</type>
      <name>getInitialConvectiveTimeSteppingType</name>
      <anchorfile>class_i_b_a_m_r_1_1_adv_diff_semi_implicit_hierarchy_integrator.html</anchorfile>
      <anchor>a2383a64024a2d6a5c263d2ece8657a48</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::CellVariable&lt; NDIM, double &gt; &gt; Q_var) const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setConvectiveOperatorType</name>
      <anchorfile>class_i_b_a_m_r_1_1_adv_diff_semi_implicit_hierarchy_integrator.html</anchorfile>
      <anchor>a8e3c3f8876eb42e6ef7a1b4289792a42</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::CellVariable&lt; NDIM, double &gt; &gt; Q_var, const std::string &amp;op_type)</arglist>
    </member>
    <member kind="function">
      <type>const std::string &amp;</type>
      <name>getConvectiveOperatorType</name>
      <anchorfile>class_i_b_a_m_r_1_1_adv_diff_semi_implicit_hierarchy_integrator.html</anchorfile>
      <anchor>a94ca2058ff9c5fd10efd5607e22ca878</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::CellVariable&lt; NDIM, double &gt; &gt; Q_var) const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setConvectiveOperatorInputDatabase</name>
      <anchorfile>class_i_b_a_m_r_1_1_adv_diff_semi_implicit_hierarchy_integrator.html</anchorfile>
      <anchor>a14cf54afef5f190afaff2f3d616bb4b6</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::CellVariable&lt; NDIM, double &gt; &gt; Q_var, SAMRAI::tbox::Pointer&lt; SAMRAI::tbox::Database &gt; input_db)</arglist>
    </member>
    <member kind="function">
      <type>SAMRAI::tbox::Pointer&lt; SAMRAI::tbox::Database &gt;</type>
      <name>getConvectiveOperatorInputDatabase</name>
      <anchorfile>class_i_b_a_m_r_1_1_adv_diff_semi_implicit_hierarchy_integrator.html</anchorfile>
      <anchor>a7d7a5c714e62f5063f180f6324c70893</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::CellVariable&lt; NDIM, double &gt; &gt; Q_var) const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setConvectiveOperator</name>
      <anchorfile>class_i_b_a_m_r_1_1_adv_diff_semi_implicit_hierarchy_integrator.html</anchorfile>
      <anchor>a580a8522f390108abcb2fc9288e69c4f</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::CellVariable&lt; NDIM, double &gt; &gt; Q_var, SAMRAI::tbox::Pointer&lt; ConvectiveOperator &gt; convective_op)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual SAMRAI::tbox::Pointer&lt; ConvectiveOperator &gt;</type>
      <name>getConvectiveOperator</name>
      <anchorfile>class_i_b_a_m_r_1_1_adv_diff_semi_implicit_hierarchy_integrator.html</anchorfile>
      <anchor>a7a269e6991b76e02f825cd004a1666ad</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::CellVariable&lt; NDIM, double &gt; &gt; Q_var)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>initializeHierarchyIntegrator</name>
      <anchorfile>class_i_b_a_m_r_1_1_adv_diff_semi_implicit_hierarchy_integrator.html</anchorfile>
      <anchor>a2cb4096a636a330a367d275c0ab60289</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::hier::PatchHierarchy&lt; NDIM &gt; &gt; hierarchy, SAMRAI::tbox::Pointer&lt; SAMRAI::mesh::GriddingAlgorithm&lt; NDIM &gt; &gt; gridding_alg)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>getNumberOfCycles</name>
      <anchorfile>class_i_b_a_m_r_1_1_adv_diff_semi_implicit_hierarchy_integrator.html</anchorfile>
      <anchor>a0e3ac85cb83ce9b6afda651335c108de</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>preprocessIntegrateHierarchy</name>
      <anchorfile>class_i_b_a_m_r_1_1_adv_diff_semi_implicit_hierarchy_integrator.html</anchorfile>
      <anchor>a7e1607b0284211c4e9f4ed483f363e93</anchor>
      <arglist>(double current_time, double new_time, int num_cycles=1)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>integrateHierarchy</name>
      <anchorfile>class_i_b_a_m_r_1_1_adv_diff_semi_implicit_hierarchy_integrator.html</anchorfile>
      <anchor>a0cea269ab13d81a65a764779ac8d2e42</anchor>
      <arglist>(double current_time, double new_time, int cycle_num=0)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>postprocessIntegrateHierarchy</name>
      <anchorfile>class_i_b_a_m_r_1_1_adv_diff_semi_implicit_hierarchy_integrator.html</anchorfile>
      <anchor>ac9b13bcce3c6af802b17c5a835323f62</anchor>
      <arglist>(double current_time, double new_time, bool skip_synchronize_new_state_data, int num_cycles=1)</arglist>
    </member>
    <member kind="function" protection="protected">
      <type>void</type>
      <name>resetHierarchyConfigurationSpecialized</name>
      <anchorfile>class_i_b_a_m_r_1_1_adv_diff_semi_implicit_hierarchy_integrator.html</anchorfile>
      <anchor>aa9391421ce827d5faa70628b0f40f2da</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::hier::BasePatchHierarchy&lt; NDIM &gt; &gt; hierarchy, int coarsest_level, int finest_level)</arglist>
    </member>
    <member kind="function" protection="protected">
      <type>void</type>
      <name>putToDatabaseSpecialized</name>
      <anchorfile>class_i_b_a_m_r_1_1_adv_diff_semi_implicit_hierarchy_integrator.html</anchorfile>
      <anchor>aca65a38a4aba2732467fc6d5b89d4745</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::tbox::Database &gt; db)</arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>TimeSteppingType</type>
      <name>d_default_convective_time_stepping_type</name>
      <anchorfile>class_i_b_a_m_r_1_1_adv_diff_semi_implicit_hierarchy_integrator.html</anchorfile>
      <anchor>ac92f7f320cbb87af122210c08837a316</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>std::string</type>
      <name>d_default_convective_op_type</name>
      <anchorfile>class_i_b_a_m_r_1_1_adv_diff_semi_implicit_hierarchy_integrator.html</anchorfile>
      <anchor>a6c1ad9780a1b2a00119a414bdd8035b0</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>std::set&lt; SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::CellVariable&lt; NDIM, double &gt; &gt; &gt;</type>
      <name>d_N_var</name>
      <anchorfile>class_i_b_a_m_r_1_1_adv_diff_semi_implicit_hierarchy_integrator.html</anchorfile>
      <anchor>ac5cb21338de1f603c73f2a08e1e2540d</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>IBAMR::AdvDiffStochasticForcing</name>
    <filename>class_i_b_a_m_r_1_1_adv_diff_stochastic_forcing.html</filename>
    <base>IBTK::CartGridFunction</base>
    <member kind="function">
      <type></type>
      <name>AdvDiffStochasticForcing</name>
      <anchorfile>class_i_b_a_m_r_1_1_adv_diff_stochastic_forcing.html</anchorfile>
      <anchor>a9c647936db6f1391da104fa80220a7d7</anchor>
      <arglist>(const std::string &amp;object_name, SAMRAI::tbox::Pointer&lt; SAMRAI::tbox::Database &gt; input_db, SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::CellVariable&lt; NDIM, double &gt; &gt; C_var, const AdvDiffSemiImplicitHierarchyIntegrator *adv_diff_solver)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~AdvDiffStochasticForcing</name>
      <anchorfile>class_i_b_a_m_r_1_1_adv_diff_stochastic_forcing.html</anchorfile>
      <anchor>a4107863c31f75505fbee99bf40ef4fe7</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>isTimeDependent</name>
      <anchorfile>class_i_b_a_m_r_1_1_adv_diff_stochastic_forcing.html</anchorfile>
      <anchor>aecd1a8e4eb5d521127643cd0b7843316</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setDataOnPatchHierarchy</name>
      <anchorfile>class_i_b_a_m_r_1_1_adv_diff_stochastic_forcing.html</anchorfile>
      <anchor>a738d91f77dcb21ad665baf87d12db11d</anchor>
      <arglist>(const int data_idx, SAMRAI::tbox::Pointer&lt; SAMRAI::hier::Variable&lt; NDIM &gt; &gt; var, SAMRAI::tbox::Pointer&lt; SAMRAI::hier::PatchHierarchy&lt; NDIM &gt; &gt; hierarchy, const double data_time, const bool initial_time=false, const int coarsest_ln=-1, const int finest_ln=-1)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setDataOnPatch</name>
      <anchorfile>class_i_b_a_m_r_1_1_adv_diff_stochastic_forcing.html</anchorfile>
      <anchor>a9804a85d5311133ac7cb835e917c0e30</anchor>
      <arglist>(const int data_idx, SAMRAI::tbox::Pointer&lt; SAMRAI::hier::Variable&lt; NDIM &gt; &gt; var, SAMRAI::tbox::Pointer&lt; SAMRAI::hier::Patch&lt; NDIM &gt; &gt; patch, const double data_time, const bool initial_time=false, SAMRAI::tbox::Pointer&lt; SAMRAI::hier::PatchLevel&lt; NDIM &gt; &gt; patch_level=SAMRAI::tbox::Pointer&lt; SAMRAI::hier::PatchLevel&lt; NDIM &gt; &gt;(NULL))</arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>std::string</type>
      <name>d_object_name</name>
      <anchorfile>class_i_b_a_m_r_1_1_adv_diff_stochastic_forcing.html</anchorfile>
      <anchor>ae640f935b97f6044f1bf486a295f763c</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>IBAMR::AdvectorExplicitPredictorPatchOps</name>
    <filename>class_i_b_a_m_r_1_1_advector_explicit_predictor_patch_ops.html</filename>
    <base>SAMRAI::tbox::Serializable</base>
    <member kind="function">
      <type></type>
      <name>AdvectorExplicitPredictorPatchOps</name>
      <anchorfile>class_i_b_a_m_r_1_1_advector_explicit_predictor_patch_ops.html</anchorfile>
      <anchor>a3fc171e433237e54643d8d85c6b3ead6</anchor>
      <arglist>(const std::string &amp;object_name, SAMRAI::tbox::Pointer&lt; SAMRAI::tbox::Database &gt; input_db, bool register_for_restart=true)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~AdvectorExplicitPredictorPatchOps</name>
      <anchorfile>class_i_b_a_m_r_1_1_advector_explicit_predictor_patch_ops.html</anchorfile>
      <anchor>af705dbef0f95265ab68534aa2dd81efa</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>const std::string &amp;</type>
      <name>getName</name>
      <anchorfile>class_i_b_a_m_r_1_1_advector_explicit_predictor_patch_ops.html</anchorfile>
      <anchor>afae2a34f72a5eb9f51953abf185c15cb</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>double</type>
      <name>computeStableDtOnPatch</name>
      <anchorfile>class_i_b_a_m_r_1_1_advector_explicit_predictor_patch_ops.html</anchorfile>
      <anchor>ac129ea3bc223f4e7d8d94f3bd790e9d6</anchor>
      <arglist>(const SAMRAI::pdat::FaceData&lt; NDIM, double &gt; &amp;u_ADV, const SAMRAI::hier::Patch&lt; NDIM &gt; &amp;patch) const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>computeAdvectiveDerivative</name>
      <anchorfile>class_i_b_a_m_r_1_1_advector_explicit_predictor_patch_ops.html</anchorfile>
      <anchor>a24f354d40171b6cc1941ef563f510a57</anchor>
      <arglist>(SAMRAI::pdat::CellData&lt; NDIM, double &gt; &amp;N, const SAMRAI::pdat::FaceData&lt; NDIM, double &gt; &amp;u_ADV, const SAMRAI::pdat::FaceData&lt; NDIM, double &gt; &amp;q_half, const SAMRAI::hier::Patch&lt; NDIM &gt; &amp;patch) const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>computeFlux</name>
      <anchorfile>class_i_b_a_m_r_1_1_advector_explicit_predictor_patch_ops.html</anchorfile>
      <anchor>a867afa66398ea4c352a264a02769e461</anchor>
      <arglist>(SAMRAI::pdat::FaceData&lt; NDIM, double &gt; &amp;flux, const SAMRAI::pdat::FaceData&lt; NDIM, double &gt; &amp;u_ADV, const SAMRAI::pdat::FaceData&lt; NDIM, double &gt; &amp;q_half, const SAMRAI::hier::Patch&lt; NDIM &gt; &amp;patch, double dt) const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>predictValue</name>
      <anchorfile>class_i_b_a_m_r_1_1_advector_explicit_predictor_patch_ops.html</anchorfile>
      <anchor>a3a703c77d69658c0ea9a820350620ee4</anchor>
      <arglist>(SAMRAI::pdat::FaceData&lt; NDIM, double &gt; &amp;q_half, const SAMRAI::pdat::FaceData&lt; NDIM, double &gt; &amp;u_ADV, const SAMRAI::pdat::CellData&lt; NDIM, double &gt; &amp;Q, const SAMRAI::hier::Patch&lt; NDIM &gt; &amp;patch, double dt) const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>predictValueWithSourceTerm</name>
      <anchorfile>class_i_b_a_m_r_1_1_advector_explicit_predictor_patch_ops.html</anchorfile>
      <anchor>a4e34ae9a617c05c24842157a1c47b964</anchor>
      <arglist>(SAMRAI::pdat::FaceData&lt; NDIM, double &gt; &amp;q_half, const SAMRAI::pdat::FaceData&lt; NDIM, double &gt; &amp;u_ADV, const SAMRAI::pdat::CellData&lt; NDIM, double &gt; &amp;Q, const SAMRAI::pdat::CellData&lt; NDIM, double &gt; &amp;F, const SAMRAI::hier::Patch&lt; NDIM &gt; &amp;patch, double dt) const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>predictNormalVelocity</name>
      <anchorfile>class_i_b_a_m_r_1_1_advector_explicit_predictor_patch_ops.html</anchorfile>
      <anchor>a92db075e4138b1fe940b319c3ee1823b</anchor>
      <arglist>(SAMRAI::pdat::FaceData&lt; NDIM, double &gt; &amp;v_half, const SAMRAI::pdat::FaceData&lt; NDIM, double &gt; &amp;u_ADV, const SAMRAI::pdat::CellData&lt; NDIM, double &gt; &amp;V, const SAMRAI::hier::Patch&lt; NDIM &gt; &amp;patch, double dt) const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>predictNormalVelocityWithSourceTerm</name>
      <anchorfile>class_i_b_a_m_r_1_1_advector_explicit_predictor_patch_ops.html</anchorfile>
      <anchor>a47d2dee2a1f708de85b8c7af3f59aa63</anchor>
      <arglist>(SAMRAI::pdat::FaceData&lt; NDIM, double &gt; &amp;v_half, const SAMRAI::pdat::FaceData&lt; NDIM, double &gt; &amp;u_ADV, const SAMRAI::pdat::CellData&lt; NDIM, double &gt; &amp;V, const SAMRAI::pdat::CellData&lt; NDIM, double &gt; &amp;F, const SAMRAI::hier::Patch&lt; NDIM &gt; &amp;patch, double dt) const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>enforceIncompressibility</name>
      <anchorfile>class_i_b_a_m_r_1_1_advector_explicit_predictor_patch_ops.html</anchorfile>
      <anchor>a858c9ea04d3da611addace65b2e65d3e</anchor>
      <arglist>(SAMRAI::pdat::FaceData&lt; NDIM, double &gt; &amp;v_half, const SAMRAI::pdat::FaceData&lt; NDIM, double &gt; &amp;u_ADV, const SAMRAI::pdat::FaceData&lt; NDIM, double &gt; &amp;grad_phi, const SAMRAI::hier::Patch&lt; NDIM &gt; &amp;patch) const </arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>getNumberCellGhosts</name>
      <anchorfile>class_i_b_a_m_r_1_1_advector_explicit_predictor_patch_ops.html</anchorfile>
      <anchor>ac07e49915a101063525f3e9f06b858ac</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>getNumberFluxGhosts</name>
      <anchorfile>class_i_b_a_m_r_1_1_advector_explicit_predictor_patch_ops.html</anchorfile>
      <anchor>a60fee32a06aa36c76e8d32c994b5fd2f</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>putToDatabase</name>
      <anchorfile>class_i_b_a_m_r_1_1_advector_explicit_predictor_patch_ops.html</anchorfile>
      <anchor>aa517fcf8700388dbac64c6eaef96fd9f</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::tbox::Database &gt; db)</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>IBAMR::AdvectorPredictorCorrectorHyperbolicPatchOps</name>
    <filename>class_i_b_a_m_r_1_1_advector_predictor_corrector_hyperbolic_patch_ops.html</filename>
    <base>HyperbolicPatchStrategy&lt; NDIM &gt;</base>
    <base>SAMRAI::tbox::Serializable</base>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>conservativeDifferenceOnPatch</name>
      <anchorfile>class_i_b_a_m_r_1_1_advector_predictor_corrector_hyperbolic_patch_ops.html</anchorfile>
      <anchor>a143d7ede1ca4d41a4e04bfc19ecba72a</anchor>
      <arglist>(SAMRAI::hier::Patch&lt; NDIM &gt; &amp;patch, double time, double dt, bool at_synchronization)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>preprocessAdvanceLevelState</name>
      <anchorfile>class_i_b_a_m_r_1_1_advector_predictor_corrector_hyperbolic_patch_ops.html</anchorfile>
      <anchor>a4dc87592761f59fd346d4b21c6492504</anchor>
      <arglist>(const SAMRAI::tbox::Pointer&lt; SAMRAI::hier::PatchLevel&lt; NDIM &gt; &gt; &amp;level, double current_time, double dt, bool first_step, bool last_step, bool regrid_advance)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>postprocessAdvanceLevelState</name>
      <anchorfile>class_i_b_a_m_r_1_1_advector_predictor_corrector_hyperbolic_patch_ops.html</anchorfile>
      <anchor>a42818a5edaeaabdd25d05968fe5ce5dc</anchor>
      <arglist>(const SAMRAI::tbox::Pointer&lt; SAMRAI::hier::PatchLevel&lt; NDIM &gt; &gt; &amp;level, double current_time, double dt, bool first_step, bool last_step, bool regrid_advance)</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>IBAMR::ConvectiveOperator</name>
    <filename>class_i_b_a_m_r_1_1_convective_operator.html</filename>
    <base>IBTK::GeneralOperator</base>
    <member kind="function" virtualness="pure">
      <type>virtual void</type>
      <name>applyConvectiveOperator</name>
      <anchorfile>class_i_b_a_m_r_1_1_convective_operator.html</anchorfile>
      <anchor>a2fef7d936e4b11450dadfbc4aa08341b</anchor>
      <arglist>(int Q_idx, int N_idx)=0</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>IBAMR::GeneralizedIBMethod</name>
    <filename>class_i_b_a_m_r_1_1_generalized_i_b_method.html</filename>
    <base>IBAMR::IBMethod</base>
    <member kind="function">
      <type></type>
      <name>GeneralizedIBMethod</name>
      <anchorfile>class_i_b_a_m_r_1_1_generalized_i_b_method.html</anchorfile>
      <anchor>abbc05e335b45a8acde97e9f435775c31</anchor>
      <arglist>(const std::string &amp;object_name, SAMRAI::tbox::Pointer&lt; SAMRAI::tbox::Database &gt; input_db, bool register_for_restart=true)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~GeneralizedIBMethod</name>
      <anchorfile>class_i_b_a_m_r_1_1_generalized_i_b_method.html</anchorfile>
      <anchor>a5e427baff82ca95291c26af1a3ef9440</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>registerIBKirchhoffRodForceGen</name>
      <anchorfile>class_i_b_a_m_r_1_1_generalized_i_b_method.html</anchorfile>
      <anchor>a6391e0e5510178a9c677866abeab1a36</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; IBKirchhoffRodForceGen &gt; ib_force_and_torque_fcn)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>registerEulerianVariables</name>
      <anchorfile>class_i_b_a_m_r_1_1_generalized_i_b_method.html</anchorfile>
      <anchor>ad8ec7bcdfb4b0b5ddcb06a847a81a6ab</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>registerEulerianCommunicationAlgorithms</name>
      <anchorfile>class_i_b_a_m_r_1_1_generalized_i_b_method.html</anchorfile>
      <anchor>a4985cb1d7efd3a3df2749cd26fbea7d4</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>preprocessIntegrateData</name>
      <anchorfile>class_i_b_a_m_r_1_1_generalized_i_b_method.html</anchorfile>
      <anchor>a8a6db2869a59241a6c1d856c9231c03c</anchor>
      <arglist>(double current_time, double new_time, int num_cycles)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>postprocessIntegrateData</name>
      <anchorfile>class_i_b_a_m_r_1_1_generalized_i_b_method.html</anchorfile>
      <anchor>a1ec24d2b9580378d377002a9c74c8ee3</anchor>
      <arglist>(double current_time, double new_time, int num_cycles)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>interpolateVelocity</name>
      <anchorfile>class_i_b_a_m_r_1_1_generalized_i_b_method.html</anchorfile>
      <anchor>ad4d7df75617ca4a9c2250c0f89ed1dc3</anchor>
      <arglist>(int u_data_idx, const std::vector&lt; SAMRAI::tbox::Pointer&lt; SAMRAI::xfer::CoarsenSchedule&lt; NDIM &gt; &gt; &gt; &amp;u_synch_scheds, const std::vector&lt; SAMRAI::tbox::Pointer&lt; SAMRAI::xfer::RefineSchedule&lt; NDIM &gt; &gt; &gt; &amp;u_ghost_fill_scheds, double data_time)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>eulerStep</name>
      <anchorfile>class_i_b_a_m_r_1_1_generalized_i_b_method.html</anchorfile>
      <anchor>ac4b221f6b8d3cef0a1ff8b0fdc0db2d9</anchor>
      <arglist>(double current_time, double new_time)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>midpointStep</name>
      <anchorfile>class_i_b_a_m_r_1_1_generalized_i_b_method.html</anchorfile>
      <anchor>a101ae8c93b783db44e0ed281b77c2215</anchor>
      <arglist>(double current_time, double new_time)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>trapezoidalStep</name>
      <anchorfile>class_i_b_a_m_r_1_1_generalized_i_b_method.html</anchorfile>
      <anchor>a4ed4ce62555ed436748e6bc8e0a552c7</anchor>
      <arglist>(double current_time, double new_time)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>computeLagrangianForce</name>
      <anchorfile>class_i_b_a_m_r_1_1_generalized_i_b_method.html</anchorfile>
      <anchor>ae82c05a6bcb5614afd9b49ca98ff725e</anchor>
      <arglist>(double data_time)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>spreadForce</name>
      <anchorfile>class_i_b_a_m_r_1_1_generalized_i_b_method.html</anchorfile>
      <anchor>aec2879c59728ba44324309bcdb6be6c9</anchor>
      <arglist>(int f_data_idx, IBTK::RobinPhysBdryPatchStrategy *f_phys_bdry_op, const std::vector&lt; SAMRAI::tbox::Pointer&lt; SAMRAI::xfer::RefineSchedule&lt; NDIM &gt; &gt; &gt; &amp;f_prolongation_scheds, double data_time)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>initializePatchHierarchy</name>
      <anchorfile>class_i_b_a_m_r_1_1_generalized_i_b_method.html</anchorfile>
      <anchor>a78a20cb6da63b19ede99e514da3b51c4</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::hier::PatchHierarchy&lt; NDIM &gt; &gt; hierarchy, SAMRAI::tbox::Pointer&lt; SAMRAI::mesh::GriddingAlgorithm&lt; NDIM &gt; &gt; gridding_alg, int u_data_idx, const std::vector&lt; SAMRAI::tbox::Pointer&lt; SAMRAI::xfer::CoarsenSchedule&lt; NDIM &gt; &gt; &gt; &amp;u_synch_scheds, const std::vector&lt; SAMRAI::tbox::Pointer&lt; SAMRAI::xfer::RefineSchedule&lt; NDIM &gt; &gt; &gt; &amp;u_ghost_fill_scheds, int integrator_step, double init_data_time, bool initial_time)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>initializeLevelData</name>
      <anchorfile>class_i_b_a_m_r_1_1_generalized_i_b_method.html</anchorfile>
      <anchor>aee56b0f43c957613f13c8fc1872a6674</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::hier::BasePatchHierarchy&lt; NDIM &gt; &gt; hierarchy, int level_number, double init_data_time, bool can_be_refined, bool initial_time, SAMRAI::tbox::Pointer&lt; SAMRAI::hier::BasePatchLevel&lt; NDIM &gt; &gt; old_level, bool allocate_data)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>putToDatabase</name>
      <anchorfile>class_i_b_a_m_r_1_1_generalized_i_b_method.html</anchorfile>
      <anchor>ae46c3894d80d3440eb1abb017829f408</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::tbox::Database &gt; db)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>IBMethod</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_method.html</anchorfile>
      <anchor>aab200ab2e593a7a9d840ccdaf9df7828</anchor>
      <arglist>(const std::string &amp;object_name, SAMRAI::tbox::Pointer&lt; SAMRAI::tbox::Database &gt; input_db, bool register_for_restart=true)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~IBMethod</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_method.html</anchorfile>
      <anchor>a2e7cbb1a8a2a78b03ecc4e514ae62500</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>registerIBLagrangianForceFunction</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_method.html</anchorfile>
      <anchor>a80be2819afc59abdee7bb79717709304</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; IBLagrangianForceStrategy &gt; ib_force_fcn)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>registerIBLagrangianSourceFunction</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_method.html</anchorfile>
      <anchor>a12fdade3f34c2bd820878628b65a3a83</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; IBLagrangianSourceStrategy &gt; ib_source_fcn)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>registerLInitStrategy</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_method.html</anchorfile>
      <anchor>ae7e2d69dca4c9e92fe9ad598033cbc69</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; IBTK::LInitStrategy &gt; l_initializer)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>freeLInitStrategy</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_method.html</anchorfile>
      <anchor>aaf75e5e82f4f1962796dd5c625257639</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>registerIBMethodPostProcessor</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_method.html</anchorfile>
      <anchor>a7a68aaa64d935433c8a0fee80d286cbf</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; IBMethodPostProcessStrategy &gt; post_processor)</arglist>
    </member>
    <member kind="function">
      <type>IBTK::LDataManager *</type>
      <name>getLDataManager</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_method.html</anchorfile>
      <anchor>a2a057348b3be36a60be7feca6f255d8f</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>SAMRAI::tbox::Pointer&lt; IBInstrumentPanel &gt;</type>
      <name>getIBInstrumentPanel</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_method.html</anchorfile>
      <anchor>a629db9f6ca60007649dce772f6ce46ae</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>registerLSiloDataWriter</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_method.html</anchorfile>
      <anchor>a42b1cd845cd4bb0ce96aa3baf6f3c683</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; IBTK::LSiloDataWriter &gt; silo_writer)</arglist>
    </member>
    <member kind="function">
      <type>const SAMRAI::hier::IntVector&lt; NDIM &gt; &amp;</type>
      <name>getMinimumGhostCellWidth</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_method.html</anchorfile>
      <anchor>aa8d10b48b12e4704e26a96beb36a286d</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setupTagBuffer</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_method.html</anchorfile>
      <anchor>a3573653e344eff869b38e2f78e979e45</anchor>
      <arglist>(SAMRAI::tbox::Array&lt; int &gt; &amp;tag_buffer, SAMRAI::tbox::Pointer&lt; SAMRAI::mesh::GriddingAlgorithm&lt; NDIM &gt; &gt; gridding_alg) const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>preprocessIntegrateData</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_method.html</anchorfile>
      <anchor>a52a9a75063ee832c621453d31cce1995</anchor>
      <arglist>(double current_time, double new_time, int num_cycles)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>postprocessIntegrateData</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_method.html</anchorfile>
      <anchor>a56b913bc4d84d06ace3f65dcd250c26c</anchor>
      <arglist>(double current_time, double new_time, int num_cycles)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>createSolverVecs</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_method.html</anchorfile>
      <anchor>a7a95029a39f946b6eded5cd296c97b51</anchor>
      <arglist>(Vec &amp;X_vec, Vec &amp;F_vec)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setupSolverVecs</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_method.html</anchorfile>
      <anchor>a8f8f0a47250913bd867d36fc1be763d1</anchor>
      <arglist>(Vec &amp;X_vec, Vec &amp;F_vec)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setUpdatedPosition</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_method.html</anchorfile>
      <anchor>a40c0494394a81e05a2b69bc46fb1e81f</anchor>
      <arglist>(Vec &amp;X_new_vec)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setLinearizedPosition</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_method.html</anchorfile>
      <anchor>a722773aeff8ee399afc997031ca60815</anchor>
      <arglist>(Vec &amp;X_vec)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>computeResidual</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_method.html</anchorfile>
      <anchor>abb195f75446d93253c0769d5c8cefb9d</anchor>
      <arglist>(Vec &amp;R_vec)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>computeLinearizedResidual</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_method.html</anchorfile>
      <anchor>acbe17e7978ffd9c8a69bfd56a4867218</anchor>
      <arglist>(Vec &amp;X_vec, Vec &amp;R_vec)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>updateFixedLEOperators</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_method.html</anchorfile>
      <anchor>a24fd286b3e46d93dcd33ef69f548d479</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>interpolateVelocity</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_method.html</anchorfile>
      <anchor>ae1ef8f1996e3412d8a93fc6713ad283d</anchor>
      <arglist>(int u_data_idx, const std::vector&lt; SAMRAI::tbox::Pointer&lt; SAMRAI::xfer::CoarsenSchedule&lt; NDIM &gt; &gt; &gt; &amp;u_synch_scheds, const std::vector&lt; SAMRAI::tbox::Pointer&lt; SAMRAI::xfer::RefineSchedule&lt; NDIM &gt; &gt; &gt; &amp;u_ghost_fill_scheds, double data_time)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>interpolateLinearizedVelocity</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_method.html</anchorfile>
      <anchor>a79e6a5bc90cdb685860d48ff018cbd78</anchor>
      <arglist>(int u_data_idx, const std::vector&lt; SAMRAI::tbox::Pointer&lt; SAMRAI::xfer::CoarsenSchedule&lt; NDIM &gt; &gt; &gt; &amp;u_synch_scheds, const std::vector&lt; SAMRAI::tbox::Pointer&lt; SAMRAI::xfer::RefineSchedule&lt; NDIM &gt; &gt; &gt; &amp;u_ghost_fill_scheds, double data_time)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>eulerStep</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_method.html</anchorfile>
      <anchor>acd7cd9797fc4a3065552cb7780ebd8dc</anchor>
      <arglist>(double current_time, double new_time)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>midpointStep</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_method.html</anchorfile>
      <anchor>a8d18c796be29076e092b3528c979bab9</anchor>
      <arglist>(double current_time, double new_time)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>trapezoidalStep</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_method.html</anchorfile>
      <anchor>ae82a3b4ffe1e325ec3ab40455ac27484</anchor>
      <arglist>(double current_time, double new_time)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>computeLagrangianForce</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_method.html</anchorfile>
      <anchor>a5f9671156c75d6c3769a939c541d7620</anchor>
      <arglist>(double data_time)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>computeLinearizedLagrangianForce</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_method.html</anchorfile>
      <anchor>ab246df397c20f8b50d2a5c8d7c1771a2</anchor>
      <arglist>(Vec &amp;X_vec, double data_time)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>spreadForce</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_method.html</anchorfile>
      <anchor>a7a47879d7390a0c73ddb1a61fbfbae34</anchor>
      <arglist>(int f_data_idx, IBTK::RobinPhysBdryPatchStrategy *f_phys_bdry_op, const std::vector&lt; SAMRAI::tbox::Pointer&lt; SAMRAI::xfer::RefineSchedule&lt; NDIM &gt; &gt; &gt; &amp;f_prolongation_scheds, double data_time)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>spreadLinearizedForce</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_method.html</anchorfile>
      <anchor>a755bdaaa0a373146f2cadaeac4f66aa9</anchor>
      <arglist>(int f_data_idx, IBTK::RobinPhysBdryPatchStrategy *f_phys_bdry_op, const std::vector&lt; SAMRAI::tbox::Pointer&lt; SAMRAI::xfer::RefineSchedule&lt; NDIM &gt; &gt; &gt; &amp;f_prolongation_scheds, double data_time)</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>hasFluidSources</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_method.html</anchorfile>
      <anchor>abb064d192e1aa401e081b3b5a4be1af6</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>computeLagrangianFluidSource</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_method.html</anchorfile>
      <anchor>a9bbe5c2dbb3bfcbaef6bf81a2f330b2c</anchor>
      <arglist>(double data_time)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>spreadFluidSource</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_method.html</anchorfile>
      <anchor>a646c17cdf4bfa09c79ce4e5160f0b6c3</anchor>
      <arglist>(int q_data_idx, const std::vector&lt; SAMRAI::tbox::Pointer&lt; SAMRAI::xfer::RefineSchedule&lt; NDIM &gt; &gt; &gt; &amp;q_prolongation_scheds, double data_time)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>interpolatePressure</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_method.html</anchorfile>
      <anchor>a19900e60018c7faa44cf17eebee9fb9f</anchor>
      <arglist>(int p_data_idx, const std::vector&lt; SAMRAI::tbox::Pointer&lt; SAMRAI::xfer::CoarsenSchedule&lt; NDIM &gt; &gt; &gt; &amp;p_synch_scheds, const std::vector&lt; SAMRAI::tbox::Pointer&lt; SAMRAI::xfer::RefineSchedule&lt; NDIM &gt; &gt; &gt; &amp;p_ghost_fill_scheds, double data_time)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>postprocessData</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_method.html</anchorfile>
      <anchor>a98ac4814cc9a96e1f66b625b7f1d98dc</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>initializePatchHierarchy</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_method.html</anchorfile>
      <anchor>adf314459e160fa2807719bea07094425</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::hier::PatchHierarchy&lt; NDIM &gt; &gt; hierarchy, SAMRAI::tbox::Pointer&lt; SAMRAI::mesh::GriddingAlgorithm&lt; NDIM &gt; &gt; gridding_alg, int u_data_idx, const std::vector&lt; SAMRAI::tbox::Pointer&lt; SAMRAI::xfer::CoarsenSchedule&lt; NDIM &gt; &gt; &gt; &amp;u_synch_scheds, const std::vector&lt; SAMRAI::tbox::Pointer&lt; SAMRAI::xfer::RefineSchedule&lt; NDIM &gt; &gt; &gt; &amp;u_ghost_fill_scheds, int integrator_step, double init_data_time, bool initial_time)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>registerLoadBalancer</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_method.html</anchorfile>
      <anchor>a933206cce0e9aed5477059f8eeee6882</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::mesh::LoadBalancer&lt; NDIM &gt; &gt; load_balancer, int workload_data_idx)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>updateWorkloadEstimates</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_method.html</anchorfile>
      <anchor>af65917771eab58e90cbd12365408db8f</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::hier::PatchHierarchy&lt; NDIM &gt; &gt; hierarchy, int workload_data_idx)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>beginDataRedistribution</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_method.html</anchorfile>
      <anchor>ad5d2b4b540e1ea1afec29095ef19d038</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::hier::PatchHierarchy&lt; NDIM &gt; &gt; hierarchy, SAMRAI::tbox::Pointer&lt; SAMRAI::mesh::GriddingAlgorithm&lt; NDIM &gt; &gt; gridding_alg)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>endDataRedistribution</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_method.html</anchorfile>
      <anchor>a6f5fca0c4895c7663af976b136481411</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::hier::PatchHierarchy&lt; NDIM &gt; &gt; hierarchy, SAMRAI::tbox::Pointer&lt; SAMRAI::mesh::GriddingAlgorithm&lt; NDIM &gt; &gt; gridding_alg)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>initializeLevelData</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_method.html</anchorfile>
      <anchor>a47183f1b49c87d0154fcf525018c5dc1</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::hier::BasePatchHierarchy&lt; NDIM &gt; &gt; hierarchy, int level_number, double init_data_time, bool can_be_refined, bool initial_time, SAMRAI::tbox::Pointer&lt; SAMRAI::hier::BasePatchLevel&lt; NDIM &gt; &gt; old_level, bool allocate_data)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>resetHierarchyConfiguration</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_method.html</anchorfile>
      <anchor>aa5449c7a669283a54fbcfb29e225971c</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::hier::BasePatchHierarchy&lt; NDIM &gt; &gt; hierarchy, int coarsest_level, int finest_level)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>applyGradientDetector</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_method.html</anchorfile>
      <anchor>a50bc8f4ac35204a8146b1650e4c35739</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::hier::BasePatchHierarchy&lt; NDIM &gt; &gt; hierarchy, int level_number, double error_data_time, int tag_index, bool initial_time, bool uses_richardson_extrapolation_too)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>putToDatabase</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_method.html</anchorfile>
      <anchor>a58e1c65f3106d63e60041342d1f26260</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::tbox::Database &gt; db)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>IBImplicitStrategy</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_implicit_strategy.html</anchorfile>
      <anchor>acc6d9b0b15130cc68744b2175f21950c</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual</type>
      <name>~IBImplicitStrategy</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_implicit_strategy.html</anchorfile>
      <anchor>a8e32352773a8afeb729b486960d4d18e</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>IBStrategy</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_strategy.html</anchorfile>
      <anchor>aef7482199cbbe0827e7653bd13d52d95</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual</type>
      <name>~IBStrategy</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_strategy.html</anchorfile>
      <anchor>a09d72676c4c86c3f4d3312e54f451d03</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>registerIBHierarchyIntegrator</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_strategy.html</anchorfile>
      <anchor>ae6818c49612138ed5a97996380acd3ea</anchor>
      <arglist>(IBHierarchyIntegrator *ib_solver)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setUseFixedLEOperators</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_strategy.html</anchorfile>
      <anchor>a7c0bace52c4d76fd7fe2ca073734a50d</anchor>
      <arglist>(bool use_fixed_coupling_ops=true)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>preprocessSolveFluidEquations</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_strategy.html</anchorfile>
      <anchor>a1f18071cd2daf62f7e51b0b872d1f2d8</anchor>
      <arglist>(double current_time, double new_time, int cycle_num)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>postprocessSolveFluidEquations</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_strategy.html</anchorfile>
      <anchor>a102efc8d6dc40f314367719889b44ed3</anchor>
      <arglist>(double current_time, double new_time, int cycle_num)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>initializeLevelData</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_strategy.html</anchorfile>
      <anchor>abeff8442684ab7303ac7a6785f516656</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::hier::BasePatchHierarchy&lt; NDIM &gt; &gt; hierarchy, int level_number, double init_data_time, bool can_be_refined, bool initial_time, SAMRAI::tbox::Pointer&lt; SAMRAI::hier::BasePatchLevel&lt; NDIM &gt; &gt; old_level, bool allocate_data)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>resetHierarchyConfiguration</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_strategy.html</anchorfile>
      <anchor>aa1802cd8fb943dde97c4bbd88642cd35</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::hier::BasePatchHierarchy&lt; NDIM &gt; &gt; hierarchy, int coarsest_level, int finest_level)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>applyGradientDetector</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_strategy.html</anchorfile>
      <anchor>a72608db38348b7ba4506aa2e2ea35482</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::hier::BasePatchHierarchy&lt; NDIM &gt; &gt; hierarchy, int level_number, double error_data_time, int tag_index, bool initial_time, bool uses_richardson_extrapolation_too)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>putToDatabase</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_strategy.html</anchorfile>
      <anchor>a285c72ab1349dd071fa2e1c0ebc3722d</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::tbox::Database &gt; db)</arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>IBHierarchyIntegrator *</type>
      <name>d_ib_solver</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_strategy.html</anchorfile>
      <anchor>aa3b7435788a1ff68e5f6f54496341011</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>bool</type>
      <name>d_use_fixed_coupling_ops</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_strategy.html</anchorfile>
      <anchor>a56c7fcd98c652f7e3f4b17e0bc59388d</anchor>
      <arglist></arglist>
    </member>
    <member kind="function" protection="protected">
      <type>void</type>
      <name>getPositionData</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_method.html</anchorfile>
      <anchor>a3918e04a644d5672f1e625e5975af6a3</anchor>
      <arglist>(std::vector&lt; SAMRAI::tbox::Pointer&lt; IBTK::LData &gt; &gt; **X_data, bool **X_needs_ghost_fill, double data_time)</arglist>
    </member>
    <member kind="function" protection="protected">
      <type>void</type>
      <name>getLinearizedPositionData</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_method.html</anchorfile>
      <anchor>ac9e1d3f8856193b21e28b1367167e5e1</anchor>
      <arglist>(std::vector&lt; SAMRAI::tbox::Pointer&lt; IBTK::LData &gt; &gt; **X_data, bool **X_needs_ghost_fill)</arglist>
    </member>
    <member kind="function" protection="protected">
      <type>void</type>
      <name>getLECouplingPositionData</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_method.html</anchorfile>
      <anchor>a20f7a43f3bc6e5f3390c186cda124517</anchor>
      <arglist>(std::vector&lt; SAMRAI::tbox::Pointer&lt; IBTK::LData &gt; &gt; **X_LE_data, bool **X_LE_needs_ghost_fill, double data_time)</arglist>
    </member>
    <member kind="function" protection="protected">
      <type>void</type>
      <name>getVelocityData</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_method.html</anchorfile>
      <anchor>a35ad74f359990098f827396f7509b381</anchor>
      <arglist>(std::vector&lt; SAMRAI::tbox::Pointer&lt; IBTK::LData &gt; &gt; **U_data, double data_time)</arglist>
    </member>
    <member kind="function" protection="protected">
      <type>void</type>
      <name>getLinearizedVelocityData</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_method.html</anchorfile>
      <anchor>ad12cc945fbdf9ea50858331d2b17b69d</anchor>
      <arglist>(std::vector&lt; SAMRAI::tbox::Pointer&lt; IBTK::LData &gt; &gt; **U_data)</arglist>
    </member>
    <member kind="function" protection="protected">
      <type>void</type>
      <name>getForceData</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_method.html</anchorfile>
      <anchor>ab4443dbb8b862b96dcaef6714bfb109d</anchor>
      <arglist>(std::vector&lt; SAMRAI::tbox::Pointer&lt; IBTK::LData &gt; &gt; **F_data, bool **F_needs_ghost_fill, double data_time)</arglist>
    </member>
    <member kind="function" protection="protected">
      <type>void</type>
      <name>getLinearizedForceData</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_method.html</anchorfile>
      <anchor>aa5e716f66635dcd1ea30e511be0b19a0</anchor>
      <arglist>(std::vector&lt; SAMRAI::tbox::Pointer&lt; IBTK::LData &gt; &gt; **F_data, bool **F_needs_ghost_fill)</arglist>
    </member>
    <member kind="function" protection="protected">
      <type>void</type>
      <name>reinitMidpointData</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_method.html</anchorfile>
      <anchor>ab050ce36e0f3284b508e02d702ccf5fd</anchor>
      <arglist>(const std::vector&lt; SAMRAI::tbox::Pointer&lt; IBTK::LData &gt; &gt; &amp;current_data, const std::vector&lt; SAMRAI::tbox::Pointer&lt; IBTK::LData &gt; &gt; &amp;new_data, const std::vector&lt; SAMRAI::tbox::Pointer&lt; IBTK::LData &gt; &gt; &amp;half_data)</arglist>
    </member>
    <member kind="function" protection="protected">
      <type>void</type>
      <name>resetAnchorPointValues</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_method.html</anchorfile>
      <anchor>a938a071fc48d5f7b32bf8f6bd4d4964d</anchor>
      <arglist>(std::vector&lt; SAMRAI::tbox::Pointer&lt; IBTK::LData &gt; &gt; U_data, int coarsest_ln, int finest_ln)</arglist>
    </member>
    <member kind="function" protection="protected">
      <type>INSHierarchyIntegrator *</type>
      <name>getINSHierarchyIntegrator</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_strategy.html</anchorfile>
      <anchor>aa36451d7a3a174554cb3199a43710b71</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function" protection="protected">
      <type>SAMRAI::tbox::Pointer&lt; SAMRAI::math::HierarchyDataOpsReal&lt; NDIM, double &gt; &gt;</type>
      <name>getVelocityHierarchyDataOps</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_strategy.html</anchorfile>
      <anchor>a7eb37654a5c2a50bb6433cd400fa4bb7</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function" protection="protected">
      <type>SAMRAI::tbox::Pointer&lt; SAMRAI::math::HierarchyDataOpsReal&lt; NDIM, double &gt; &gt;</type>
      <name>getPressureHierarchyDataOps</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_strategy.html</anchorfile>
      <anchor>a3c10c980dd56fe9c345e877647f33cd3</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function" protection="protected">
      <type>SAMRAI::tbox::Pointer&lt; IBTK::HierarchyMathOps &gt;</type>
      <name>getHierarchyMathOps</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_strategy.html</anchorfile>
      <anchor>a019ef4427f51459508cb603f3d0c7c84</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function" protection="protected">
      <type>void</type>
      <name>registerVariable</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_strategy.html</anchorfile>
      <anchor>a9eb305eda7681c70e82dd18aba06535f</anchor>
      <arglist>(int &amp;current_idx, int &amp;new_idx, int &amp;scratch_idx, SAMRAI::tbox::Pointer&lt; SAMRAI::hier::Variable&lt; NDIM &gt; &gt; variable, const SAMRAI::hier::IntVector&lt; NDIM &gt; &amp;scratch_ghosts=SAMRAI::hier::IntVector&lt; NDIM &gt;(0), const std::string &amp;coarsen_name=&quot;NO_COARSEN&quot;, const std::string &amp;refine_name=&quot;NO_REFINE&quot;, SAMRAI::tbox::Pointer&lt; IBTK::CartGridFunction &gt; init_fcn=SAMRAI::tbox::Pointer&lt; IBTK::CartGridFunction &gt;(NULL))</arglist>
    </member>
    <member kind="function" protection="protected">
      <type>void</type>
      <name>registerVariable</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_strategy.html</anchorfile>
      <anchor>ac6b30572b0b5dacbbe72b4f6a2531512</anchor>
      <arglist>(int &amp;idx, SAMRAI::tbox::Pointer&lt; SAMRAI::hier::Variable&lt; NDIM &gt; &gt; variable, const SAMRAI::hier::IntVector&lt; NDIM &gt; &amp;ghosts=SAMRAI::hier::IntVector&lt; NDIM &gt;(0), SAMRAI::tbox::Pointer&lt; SAMRAI::hier::VariableContext &gt; ctx=SAMRAI::tbox::Pointer&lt; SAMRAI::hier::VariableContext &gt;(NULL))</arglist>
    </member>
    <member kind="function" protection="protected">
      <type>void</type>
      <name>registerGhostfillRefineAlgorithm</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_strategy.html</anchorfile>
      <anchor>a5e30d3ee18e368de135a89677edfa279</anchor>
      <arglist>(const std::string &amp;name, SAMRAI::tbox::Pointer&lt; SAMRAI::xfer::RefineAlgorithm&lt; NDIM &gt; &gt; ghostfill_alg, SAMRAI::xfer::RefinePatchStrategy&lt; NDIM &gt; *ghostfill_patch_strategy=NULL)</arglist>
    </member>
    <member kind="function" protection="protected">
      <type>void</type>
      <name>registerProlongRefineAlgorithm</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_strategy.html</anchorfile>
      <anchor>aa42fcb8aa8327dfc39bde69bd7462a68</anchor>
      <arglist>(const std::string &amp;name, SAMRAI::tbox::Pointer&lt; SAMRAI::xfer::RefineAlgorithm&lt; NDIM &gt; &gt; prolong_alg, SAMRAI::xfer::RefinePatchStrategy&lt; NDIM &gt; *prolong_patch_strategy=NULL)</arglist>
    </member>
    <member kind="function" protection="protected">
      <type>void</type>
      <name>registerCoarsenAlgorithm</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_strategy.html</anchorfile>
      <anchor>a5ef8e340438b75e855a3c81018251286</anchor>
      <arglist>(const std::string &amp;name, SAMRAI::tbox::Pointer&lt; SAMRAI::xfer::CoarsenAlgorithm&lt; NDIM &gt; &gt; coarsen_alg, SAMRAI::xfer::CoarsenPatchStrategy&lt; NDIM &gt; *coarsen_patch_strategy=NULL)</arglist>
    </member>
    <member kind="function" protection="protected">
      <type>SAMRAI::tbox::Pointer&lt; SAMRAI::xfer::RefineAlgorithm&lt; NDIM &gt; &gt;</type>
      <name>getGhostfillRefineAlgorithm</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_strategy.html</anchorfile>
      <anchor>a4f351325a0e6f2a23daa862a65273199</anchor>
      <arglist>(const std::string &amp;name) const </arglist>
    </member>
    <member kind="function" protection="protected">
      <type>SAMRAI::tbox::Pointer&lt; SAMRAI::xfer::RefineAlgorithm&lt; NDIM &gt; &gt;</type>
      <name>getProlongRefineAlgorithm</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_strategy.html</anchorfile>
      <anchor>a72b0a99686049e6b97a86d63c73df496</anchor>
      <arglist>(const std::string &amp;name) const </arglist>
    </member>
    <member kind="function" protection="protected">
      <type>SAMRAI::tbox::Pointer&lt; SAMRAI::xfer::CoarsenAlgorithm&lt; NDIM &gt; &gt;</type>
      <name>getCoarsenAlgorithm</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_strategy.html</anchorfile>
      <anchor>a5d96b6e42463174d23e98b6f0fd66367</anchor>
      <arglist>(const std::string &amp;name) const </arglist>
    </member>
    <member kind="function" protection="protected">
      <type>const std::vector&lt; SAMRAI::tbox::Pointer&lt; SAMRAI::xfer::RefineSchedule&lt; NDIM &gt; &gt; &gt; &amp;</type>
      <name>getGhostfillRefineSchedules</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_strategy.html</anchorfile>
      <anchor>ab2b6abeaeac716e6ef9767a52aefab76</anchor>
      <arglist>(const std::string &amp;name) const </arglist>
    </member>
    <member kind="function" protection="protected">
      <type>const std::vector&lt; SAMRAI::tbox::Pointer&lt; SAMRAI::xfer::RefineSchedule&lt; NDIM &gt; &gt; &gt; &amp;</type>
      <name>getProlongRefineSchedules</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_strategy.html</anchorfile>
      <anchor>a3b9dd5359d7809b8217a5f3ab99fde68</anchor>
      <arglist>(const std::string &amp;name) const </arglist>
    </member>
    <member kind="function" protection="protected">
      <type>const std::vector&lt; SAMRAI::tbox::Pointer&lt; SAMRAI::xfer::CoarsenSchedule&lt; NDIM &gt; &gt; &gt; &amp;</type>
      <name>getCoarsenSchedules</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_strategy.html</anchorfile>
      <anchor>a39de92fb247539e7ddd3b621521467e7</anchor>
      <arglist>(const std::string &amp;name) const </arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>IBAMR::IBAnchorPointSpec</name>
    <filename>class_i_b_a_m_r_1_1_i_b_anchor_point_spec.html</filename>
    <base>IBTK::Streamable</base>
    <member kind="function">
      <type></type>
      <name>IBAnchorPointSpec</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_anchor_point_spec.html</anchorfile>
      <anchor>a85ef5fa605ad79ecf9b34756638963e4</anchor>
      <arglist>(int node_idx=-1)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~IBAnchorPointSpec</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_anchor_point_spec.html</anchorfile>
      <anchor>aa75d7127c45bbc54f11552d5aa5d3653</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>const int &amp;</type>
      <name>getNodeIndex</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_anchor_point_spec.html</anchorfile>
      <anchor>aa940a5df3abc8374b5509b69e9eb44be</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>int &amp;</type>
      <name>getNodeIndex</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_anchor_point_spec.html</anchorfile>
      <anchor>a282bd8fb3382a6b58d33aea001deedea</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>getStreamableClassID</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_anchor_point_spec.html</anchorfile>
      <anchor>ac7c6a125cf5c8ea76613a4be7084aafb</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>size_t</type>
      <name>getDataStreamSize</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_anchor_point_spec.html</anchorfile>
      <anchor>ac33e99ff3f9e5f1030f467737d1daa50</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>packStream</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_anchor_point_spec.html</anchorfile>
      <anchor>a1b04b5b8fdb6271685ec15afb146b130</anchor>
      <arglist>(SAMRAI::tbox::AbstractStream &amp;stream)</arglist>
    </member>
    <member kind="function" static="yes">
      <type>static void</type>
      <name>registerWithStreamableManager</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_anchor_point_spec.html</anchorfile>
      <anchor>accfda36202da1ba08123558c0c4ff0c7</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" static="yes">
      <type>static bool</type>
      <name>getIsRegisteredWithStreamableManager</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_anchor_point_spec.html</anchorfile>
      <anchor>a2ace9c93d379007f2765dbc12115ae5d</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="variable" static="yes">
      <type>static int</type>
      <name>STREAMABLE_CLASS_ID</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_anchor_point_spec.html</anchorfile>
      <anchor>a4b3ca38906344c0ffcaf00f19064bc36</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>IBAMR::IBBeamForceSpec</name>
    <filename>class_i_b_a_m_r_1_1_i_b_beam_force_spec.html</filename>
    <base>IBTK::Streamable</base>
    <member kind="typedef">
      <type>std::pair&lt; int, int &gt;</type>
      <name>NeighborIdxs</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_beam_force_spec.html</anchorfile>
      <anchor>a16dba2c366880ada662515a2e02da607</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>IBBeamForceSpec</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_beam_force_spec.html</anchorfile>
      <anchor>ad1c2b8bf45fc6c8a5425939e7b6524ae</anchor>
      <arglist>(unsigned int num_beams=0)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>IBBeamForceSpec</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_beam_force_spec.html</anchorfile>
      <anchor>a979f4fb444f99cfce54fd0b82df15814</anchor>
      <arglist>(int master_idx, const std::vector&lt; NeighborIdxs &gt; &amp;neighbor_idxs, const std::vector&lt; double &gt; &amp;bend_rigidities, const std::vector&lt; IBTK::Vector &gt; &amp;mesh_dependent_curvatures)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~IBBeamForceSpec</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_beam_force_spec.html</anchorfile>
      <anchor>a41a0b29477062712995413ed2f1d7207</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>unsigned int</type>
      <name>getNumberOfBeams</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_beam_force_spec.html</anchorfile>
      <anchor>ab83502ca1a6087acbca43dfc75b3415c</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>const int &amp;</type>
      <name>getMasterNodeIndex</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_beam_force_spec.html</anchorfile>
      <anchor>af22af2c07a7987c6facb4d5522a834be</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>int &amp;</type>
      <name>getMasterNodeIndex</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_beam_force_spec.html</anchorfile>
      <anchor>ac1fcf122ecff666fe5bd03b9a4ae9c28</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>const std::vector&lt; NeighborIdxs &gt; &amp;</type>
      <name>getNeighborNodeIndices</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_beam_force_spec.html</anchorfile>
      <anchor>a07e76e17387190ef0fbcf490769549c3</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>std::vector&lt; NeighborIdxs &gt; &amp;</type>
      <name>getNeighborNodeIndices</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_beam_force_spec.html</anchorfile>
      <anchor>a63a34a75f65c8202f0f574bc0f9f8244</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>const std::vector&lt; double &gt; &amp;</type>
      <name>getBendingRigidities</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_beam_force_spec.html</anchorfile>
      <anchor>ac7bd19c036e4a0fb2af7cc0780468df3</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>std::vector&lt; double &gt; &amp;</type>
      <name>getBendingRigidities</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_beam_force_spec.html</anchorfile>
      <anchor>aebc008b69f54721d52695b3b2df793a7</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>const std::vector&lt; IBTK::Vector &gt; &amp;</type>
      <name>getMeshDependentCurvatures</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_beam_force_spec.html</anchorfile>
      <anchor>a02b78a9672d0fa647a02b2b42b1ad375</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>std::vector&lt; IBTK::Vector &gt; &amp;</type>
      <name>getMeshDependentCurvatures</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_beam_force_spec.html</anchorfile>
      <anchor>a55e7e7a9495cbdbebbd2538987d34d61</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>getStreamableClassID</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_beam_force_spec.html</anchorfile>
      <anchor>ac0715e93d0c961978900ab6a483df88f</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>size_t</type>
      <name>getDataStreamSize</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_beam_force_spec.html</anchorfile>
      <anchor>aeec173447c7ceb7c4fb8e6906a11164d</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>packStream</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_beam_force_spec.html</anchorfile>
      <anchor>aba9dc85b401fcdef7b84af4ba1dd8dec</anchor>
      <arglist>(SAMRAI::tbox::AbstractStream &amp;stream)</arglist>
    </member>
    <member kind="function" static="yes">
      <type>static void</type>
      <name>registerWithStreamableManager</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_beam_force_spec.html</anchorfile>
      <anchor>afdb8f87dd2ee4637e13cf570387a2fff</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" static="yes">
      <type>static bool</type>
      <name>getIsRegisteredWithStreamableManager</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_beam_force_spec.html</anchorfile>
      <anchor>ad76ae0320c080659f518aa7606ebf011</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="variable" static="yes">
      <type>static int</type>
      <name>STREAMABLE_CLASS_ID</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_beam_force_spec.html</anchorfile>
      <anchor>aab01a835427b9e03f8db0aabd6987427</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>IBAMR::IBExplicitHierarchyIntegrator</name>
    <filename>class_i_b_a_m_r_1_1_i_b_explicit_hierarchy_integrator.html</filename>
    <base>IBAMR::IBHierarchyIntegrator</base>
    <member kind="function">
      <type></type>
      <name>IBExplicitHierarchyIntegrator</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_explicit_hierarchy_integrator.html</anchorfile>
      <anchor>a0e009f63e9382331f5e137bc24518a0e</anchor>
      <arglist>(const std::string &amp;object_name, SAMRAI::tbox::Pointer&lt; SAMRAI::tbox::Database &gt; input_db, SAMRAI::tbox::Pointer&lt; IBStrategy &gt; ib_method_ops, SAMRAI::tbox::Pointer&lt; INSHierarchyIntegrator &gt; ins_hier_integrator, bool register_for_restart=true)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~IBExplicitHierarchyIntegrator</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_explicit_hierarchy_integrator.html</anchorfile>
      <anchor>add3a4bfd2ad5e1b56fff5e84f4430860</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>preprocessIntegrateHierarchy</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_explicit_hierarchy_integrator.html</anchorfile>
      <anchor>ab703f58ac4ec7d089508ca5ae5f053ff</anchor>
      <arglist>(double current_time, double new_time, int num_cycles=1)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>integrateHierarchy</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_explicit_hierarchy_integrator.html</anchorfile>
      <anchor>a3be22787837595e9884bb7766ff8f260</anchor>
      <arglist>(double current_time, double new_time, int cycle_num=0)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>postprocessIntegrateHierarchy</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_explicit_hierarchy_integrator.html</anchorfile>
      <anchor>a62bf65e268882dbfb334f970a6442d19</anchor>
      <arglist>(double current_time, double new_time, bool skip_synchronize_new_state_data, int num_cycles=1)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>initializeHierarchyIntegrator</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_explicit_hierarchy_integrator.html</anchorfile>
      <anchor>a25e401819d493c5cbfbcde9c719d9f05</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::hier::PatchHierarchy&lt; NDIM &gt; &gt; hierarchy, SAMRAI::tbox::Pointer&lt; SAMRAI::mesh::GriddingAlgorithm&lt; NDIM &gt; &gt; gridding_alg)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~IBHierarchyIntegrator</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_hierarchy_integrator.html</anchorfile>
      <anchor>a046dbd07547d62f8591570d84a8720a5</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>SAMRAI::tbox::Pointer&lt; IBStrategy &gt;</type>
      <name>getIBStrategy</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_hierarchy_integrator.html</anchorfile>
      <anchor>a3fa92794c0eb0a0be4aefc2b27f31cee</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>registerBodyForceFunction</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_hierarchy_integrator.html</anchorfile>
      <anchor>a89cda645ab27a9a4e02cb17b9042a022</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; IBTK::CartGridFunction &gt; F_fcn)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>registerLoadBalancer</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_hierarchy_integrator.html</anchorfile>
      <anchor>aa4c611fb549acc746272e9b929c83852</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::mesh::LoadBalancer&lt; NDIM &gt; &gt; load_balancer)</arglist>
    </member>
    <member kind="function">
      <type>SAMRAI::tbox::Pointer&lt; SAMRAI::hier::Variable&lt; NDIM &gt; &gt;</type>
      <name>getVelocityVariable</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_hierarchy_integrator.html</anchorfile>
      <anchor>acd6a0977bbe750bd9ec4a6dbdff05d43</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>SAMRAI::tbox::Pointer&lt; SAMRAI::hier::Variable&lt; NDIM &gt; &gt;</type>
      <name>getPressureVariable</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_hierarchy_integrator.html</anchorfile>
      <anchor>a9cb8f415fdcbd4a21e9d235c8060f855</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>SAMRAI::tbox::Pointer&lt; SAMRAI::hier::Variable&lt; NDIM &gt; &gt;</type>
      <name>getBodyForceVariable</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_hierarchy_integrator.html</anchorfile>
      <anchor>a0c888c97e0e85c2ff43b976722609596</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>SAMRAI::tbox::Pointer&lt; SAMRAI::hier::Variable&lt; NDIM &gt; &gt;</type>
      <name>getFluidSourceVariable</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_hierarchy_integrator.html</anchorfile>
      <anchor>a32bfb45380a55ccacde20e341d3e5ee5</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>initializeHierarchyIntegrator</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_hierarchy_integrator.html</anchorfile>
      <anchor>afb457157586f1ed77fecd357dcf20d8b</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::hier::PatchHierarchy&lt; NDIM &gt; &gt; hierarchy, SAMRAI::tbox::Pointer&lt; SAMRAI::mesh::GriddingAlgorithm&lt; NDIM &gt; &gt; gridding_alg)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>initializePatchHierarchy</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_hierarchy_integrator.html</anchorfile>
      <anchor>ad4569b55cf00263ebb79f7a79792b656</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::hier::PatchHierarchy&lt; NDIM &gt; &gt; hierarchy, SAMRAI::tbox::Pointer&lt; SAMRAI::mesh::GriddingAlgorithm&lt; NDIM &gt; &gt; gridding_alg)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>regridHierarchy</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_hierarchy_integrator.html</anchorfile>
      <anchor>a79ab4eddf2f3757fc5fb171f0dd9e264</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" protection="protected">
      <type>void</type>
      <name>putToDatabaseSpecialized</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_explicit_hierarchy_integrator.html</anchorfile>
      <anchor>ac1273f0917b438293d7e1a41388fbea0</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::tbox::Database &gt; db)</arglist>
    </member>
    <member kind="function" protection="protected">
      <type></type>
      <name>IBHierarchyIntegrator</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_hierarchy_integrator.html</anchorfile>
      <anchor>a715f808b9523ccbf846c4ca9bc399782</anchor>
      <arglist>(const std::string &amp;object_name, SAMRAI::tbox::Pointer&lt; SAMRAI::tbox::Database &gt; input_db, SAMRAI::tbox::Pointer&lt; IBStrategy &gt; ib_method_ops, SAMRAI::tbox::Pointer&lt; INSHierarchyIntegrator &gt; ins_hier_integrator, bool register_for_restart=true)</arglist>
    </member>
    <member kind="function" protection="protected">
      <type>bool</type>
      <name>atRegridPointSpecialized</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_hierarchy_integrator.html</anchorfile>
      <anchor>a1e32170f574fb296b1543182583dcce9</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function" protection="protected">
      <type>void</type>
      <name>initializeLevelDataSpecialized</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_hierarchy_integrator.html</anchorfile>
      <anchor>add9534679340a51533d7db7c895894a9</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::hier::BasePatchHierarchy&lt; NDIM &gt; &gt; hierarchy, int level_number, double init_data_time, bool can_be_refined, bool initial_time, SAMRAI::tbox::Pointer&lt; SAMRAI::hier::BasePatchLevel&lt; NDIM &gt; &gt; old_level, bool allocate_data)</arglist>
    </member>
    <member kind="function" protection="protected">
      <type>void</type>
      <name>resetHierarchyConfigurationSpecialized</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_hierarchy_integrator.html</anchorfile>
      <anchor>a0ba04f93db924d6f1580acb0f054866b</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::hier::BasePatchHierarchy&lt; NDIM &gt; &gt; hierarchy, int coarsest_level, int finest_level)</arglist>
    </member>
    <member kind="function" protection="protected">
      <type>void</type>
      <name>applyGradientDetectorSpecialized</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_hierarchy_integrator.html</anchorfile>
      <anchor>af3c702921b0c37be050ca2f8f6d84b28</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::hier::BasePatchHierarchy&lt; NDIM &gt; &gt; hierarchy, int level_number, double error_data_time, int tag_index, bool initial_time, bool uses_richardson_extrapolation_too)</arglist>
    </member>
    <member kind="function" protection="protected">
      <type>void</type>
      <name>putToDatabaseSpecialized</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_hierarchy_integrator.html</anchorfile>
      <anchor>ac1d745eddb7f546afcebc9b07384a83f</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::tbox::Database &gt; db)</arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>TimeSteppingType</type>
      <name>d_time_stepping_type</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_hierarchy_integrator.html</anchorfile>
      <anchor>af6e3966b28f211928949a659e7b77572</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>bool</type>
      <name>d_error_on_dt_change</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_hierarchy_integrator.html</anchorfile>
      <anchor>a0ae88f36249acfb3058bfe9b2c0d5c5f</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>IBAMR::IBFECentroidPostProcessor</name>
    <filename>class_i_b_a_m_r_1_1_i_b_f_e_centroid_post_processor.html</filename>
    <base>IBAMR::IBFEPostProcessor</base>
    <member kind="function">
      <type></type>
      <name>IBFECentroidPostProcessor</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_f_e_centroid_post_processor.html</anchorfile>
      <anchor>ac34986df67a6d7f3ceff670c3e3f3d44</anchor>
      <arglist>(const std::string &amp;name, IBTK::FEDataManager *fe_data_manager)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~IBFECentroidPostProcessor</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_f_e_centroid_post_processor.html</anchorfile>
      <anchor>a3ed987508bb73857ba0b9a2b84db1d5f</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>registerScalarVariable</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_f_e_centroid_post_processor.html</anchorfile>
      <anchor>a30db6e6ee19a97477d907ddc1246e61d</anchor>
      <arglist>(const std::string &amp;var_name, libMeshEnums::FEFamily var_fe_family, libMeshEnums::Order var_fe_order, IBTK::ScalarMeshFcnPtr var_fcn, std::vector&lt; unsigned int &gt; var_fcn_systems=std::vector&lt; unsigned int &gt;(), void *var_fcn_ctx=NULL)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>registerVectorVariable</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_f_e_centroid_post_processor.html</anchorfile>
      <anchor>ac8c2229d73ecf502f97895a51d1c5360</anchor>
      <arglist>(const std::string &amp;var_name, libMeshEnums::FEFamily var_fe_family, libMeshEnums::Order var_fe_order, IBTK::VectorMeshFcnPtr var_fcn, std::vector&lt; unsigned int &gt; var_fcn_systems=std::vector&lt; unsigned int &gt;(), void *var_fcn_ctx=NULL, unsigned int var_dim=NDIM)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>registerTensorVariable</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_f_e_centroid_post_processor.html</anchorfile>
      <anchor>a76bf816be2e404758c4ec41ea7454801</anchor>
      <arglist>(const std::string &amp;var_name, libMeshEnums::FEFamily var_fe_family, libMeshEnums::Order var_fe_order, IBTK::TensorMeshFcnPtr var_fcn, std::vector&lt; unsigned int &gt; var_fcn_systems=std::vector&lt; unsigned int &gt;(), void *var_fcn_ctx=NULL, unsigned int var_dim=NDIM)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>reconstructVariables</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_f_e_centroid_post_processor.html</anchorfile>
      <anchor>ad3912686c68254ba1ffcba041020ba28</anchor>
      <arglist>(double data_time)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>IBFEPostProcessor</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_f_e_post_processor.html</anchorfile>
      <anchor>a0f4748051d1c23defcde01303e0ab15c</anchor>
      <arglist>(const std::string &amp;name, IBTK::FEDataManager *fe_data_manager)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual</type>
      <name>~IBFEPostProcessor</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_f_e_post_processor.html</anchorfile>
      <anchor>aab0050f2261e63e404dbce50543e8942</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>registerInterpolatedScalarEulerianVariable</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_f_e_post_processor.html</anchorfile>
      <anchor>a2dbe40e9087d5f843cf463d3bdb7a15c</anchor>
      <arglist>(const std::string &amp;var_name, libMeshEnums::FEFamily var_fe_family, libMeshEnums::Order var_fe_order, SAMRAI::tbox::Pointer&lt; SAMRAI::hier::Variable&lt; NDIM &gt; &gt; var, SAMRAI::tbox::Pointer&lt; SAMRAI::hier::VariableContext &gt; ctx)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>registerInterpolatedScalarEulerianVariable</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_f_e_post_processor.html</anchorfile>
      <anchor>a4492741d8d5e42c90ef614c33ab489dc</anchor>
      <arglist>(const std::string &amp;var_name, libMeshEnums::FEFamily var_fe_family, libMeshEnums::Order var_fe_order, SAMRAI::tbox::Pointer&lt; SAMRAI::hier::Variable&lt; NDIM &gt; &gt; var, SAMRAI::tbox::Pointer&lt; SAMRAI::hier::VariableContext &gt; ctx, const IBTK::FEDataManager::InterpSpec &amp;interp_spec)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>initializeFEData</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_f_e_post_processor.html</anchorfile>
      <anchor>a2356dce6c5e46014ca9151c858177a55</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>postProcessData</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_f_e_post_processor.html</anchorfile>
      <anchor>a379d2d1a553d2c89211970f4d600f775</anchor>
      <arglist>(double data_time)</arglist>
    </member>
    <member kind="function" static="yes">
      <type>static void</type>
      <name>FF_fcn</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_f_e_post_processor.html</anchorfile>
      <anchor>ab994a8a6263c12471f2101a82d90bc82</anchor>
      <arglist>(libMesh::TensorValue&lt; double &gt; &amp;FF_out, const libMesh::TensorValue&lt; double &gt; &amp;FF_in, const libMesh::Point &amp;, const libMesh::Point &amp;, libMesh::Elem *, const std::vector&lt; libMesh::NumericVector&lt; double &gt; * &gt; &amp;, double, void *)</arglist>
    </member>
    <member kind="function" static="yes">
      <type>static void</type>
      <name>EE_fcn</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_f_e_post_processor.html</anchorfile>
      <anchor>a8e29e4a7dbcdbb955c3860aea1b79eb2</anchor>
      <arglist>(libMesh::TensorValue&lt; double &gt; &amp;EE, const libMesh::TensorValue&lt; double &gt; &amp;FF, const libMesh::Point &amp;, const libMesh::Point &amp;, libMesh::Elem *, const std::vector&lt; libMesh::NumericVector&lt; double &gt; * &gt; &amp;, double, void *)</arglist>
    </member>
    <member kind="function" static="yes">
      <type>static void</type>
      <name>cauchy_stress_from_PK1_stress_fcn</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_f_e_post_processor.html</anchorfile>
      <anchor>a0afc94ad949693e74edbb9be45f7d428</anchor>
      <arglist>(libMesh::TensorValue&lt; double &gt; &amp;sigma, const libMesh::TensorValue&lt; double &gt; &amp;FF, const libMesh::Point &amp;X, const libMesh::Point &amp;s, libMesh::Elem *elem, const std::vector&lt; libMesh::NumericVector&lt; double &gt; * &gt; &amp;system_data, double data_time, void *ctx)</arglist>
    </member>
    <member kind="function" static="yes">
      <type>static void</type>
      <name>deformed_material_axis_fcn</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_f_e_post_processor.html</anchorfile>
      <anchor>ab9f4d13cd2990d15a66dd8f40c880f52</anchor>
      <arglist>(libMesh::VectorValue&lt; double &gt; &amp;f, const libMesh::TensorValue&lt; double &gt; &amp;FF, const libMesh::Point &amp;, const libMesh::Point &amp;, libMesh::Elem *elem, const std::vector&lt; libMesh::NumericVector&lt; double &gt; * &gt; &amp;system_data, double, void *ctx)</arglist>
    </member>
    <member kind="function" static="yes">
      <type>static void</type>
      <name>deformed_normalized_material_axis_fcn</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_f_e_post_processor.html</anchorfile>
      <anchor>a8d45b889f23bfede779f5ea89fefe41e</anchor>
      <arglist>(libMesh::VectorValue&lt; double &gt; &amp;f, const libMesh::TensorValue&lt; double &gt; &amp;FF, const libMesh::Point &amp;, const libMesh::Point &amp;, libMesh::Elem *elem, const std::vector&lt; libMesh::NumericVector&lt; double &gt; * &gt; &amp;system_data, double, void *ctx)</arglist>
    </member>
    <member kind="function" static="yes">
      <type>static void</type>
      <name>material_axis_stretch_fcn</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_f_e_post_processor.html</anchorfile>
      <anchor>abef6156d8dccbffed11e248b00b1bf0c</anchor>
      <arglist>(double &amp;lambda, const libMesh::TensorValue&lt; double &gt; &amp;FF, const libMesh::Point &amp;, const libMesh::Point &amp;, libMesh::Elem *elem, const std::vector&lt; libMesh::NumericVector&lt; double &gt; * &gt; &amp;system_data, double, void *ctx)</arglist>
    </member>
    <member kind="function" protection="protected" virtualness="virtual">
      <type>virtual void</type>
      <name>interpolateVariables</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_f_e_post_processor.html</anchorfile>
      <anchor>a3d9e3d89f208a3280428709c2cc08bb7</anchor>
      <arglist>(double data_time)</arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>const std::string</type>
      <name>d_name</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_f_e_post_processor.html</anchorfile>
      <anchor>a56c844fad1e9e9bfec0409afacacd418</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>libMesh::MeshBase *</type>
      <name>d_mesh</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_f_e_post_processor.html</anchorfile>
      <anchor>a185d12c9d948fa74589a556bb5b71fb3</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>std::vector&lt; libMesh::System * &gt;</type>
      <name>d_scalar_var_systems</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_f_e_post_processor.html</anchorfile>
      <anchor>a8258793f4317c7da38602a82117ec133</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>std::vector&lt; libMesh::System * &gt;</type>
      <name>d_vector_var_systems</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_f_e_post_processor.html</anchorfile>
      <anchor>afe45f871a6f47b0fdb619b028d27e717</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>std::vector&lt; libMesh::System * &gt;</type>
      <name>d_tensor_var_systems</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_f_e_post_processor.html</anchorfile>
      <anchor>ae399fd9c76ed72af006faaabab21a135</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>std::vector&lt; libMesh::System * &gt;</type>
      <name>d_scalar_interp_var_systems</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_f_e_post_processor.html</anchorfile>
      <anchor>acb143873916e5ee50e9205f41ab09ef2</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>std::vector&lt; libMesh::System * &gt;</type>
      <name>d_var_systems</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_f_e_post_processor.html</anchorfile>
      <anchor>a2f923b373bde32b3a6152c67c3022011</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>std::set&lt; unsigned int &gt;</type>
      <name>d_var_fcn_systems</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_f_e_post_processor.html</anchorfile>
      <anchor>a06805b382b6fdcf85c015722c44d76c9</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>IBAMR::IBFEMethod</name>
    <filename>class_i_b_a_m_r_1_1_i_b_f_e_method.html</filename>
    <base>IBAMR::IBStrategy</base>
    <class kind="struct">IBAMR::IBFEMethod::ConstrainedVelocityFcnData</class>
    <class kind="struct">IBAMR::IBFEMethod::CoordinateMappingFcnData</class>
    <class kind="struct">IBAMR::IBFEMethod::LagBodyForceFcnData</class>
    <class kind="struct">IBAMR::IBFEMethod::LagSurfaceForceFcnData</class>
    <class kind="struct">IBAMR::IBFEMethod::LagSurfacePressureFcnData</class>
    <class kind="struct">IBAMR::IBFEMethod::PK1StressFcnData</class>
    <member kind="typedef">
      <type>void(*</type>
      <name>ConstrainedVelocityFcnPtr</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_f_e_method.html</anchorfile>
      <anchor>a8d5d06f0f170003d8e048174f4396a50</anchor>
      <arglist>)(libMesh::NumericVector&lt; double &gt; &amp;U_b, libMesh::NumericVector&lt; double &gt; &amp;U, libMesh::NumericVector&lt; double &gt; &amp;X, libMesh::EquationSystems *equation_systems, double data_time, void *ctx)</arglist>
    </member>
    <member kind="typedef">
      <type>void(*</type>
      <name>CoordinateMappingFcnPtr</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_f_e_method.html</anchorfile>
      <anchor>a08824abcee1f911158f0721a1d7514cf</anchor>
      <arglist>)(libMesh::Point &amp;X, const libMesh::Point &amp;s, void *ctx)</arglist>
    </member>
    <member kind="typedef">
      <type>IBTK::TensorMeshFcnPtr</type>
      <name>PK1StressFcnPtr</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_f_e_method.html</anchorfile>
      <anchor>a3e238b84126b64acbac30dca07312eee</anchor>
      <arglist></arglist>
    </member>
    <member kind="typedef">
      <type>IBTK::VectorMeshFcnPtr</type>
      <name>LagBodyForceFcnPtr</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_f_e_method.html</anchorfile>
      <anchor>a3f9cbf889dbf4ec311d4d7d7e062830d</anchor>
      <arglist></arglist>
    </member>
    <member kind="typedef">
      <type>IBTK::ScalarSurfaceFcnPtr</type>
      <name>LagSurfacePressureFcnPtr</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_f_e_method.html</anchorfile>
      <anchor>adeb5cb97129bc2debd0741a858b41f73</anchor>
      <arglist></arglist>
    </member>
    <member kind="typedef">
      <type>IBTK::VectorSurfaceFcnPtr</type>
      <name>LagSurfaceForceFcnPtr</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_f_e_method.html</anchorfile>
      <anchor>a441dedcd4e9b6f94ddc36d6b5afd7c4a</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>IBFEMethod</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_f_e_method.html</anchorfile>
      <anchor>a36bcb0c6217c56810875e8c5bd85a464</anchor>
      <arglist>(const std::string &amp;object_name, SAMRAI::tbox::Pointer&lt; SAMRAI::tbox::Database &gt; input_db, libMesh::Mesh *mesh, int max_level_number, bool register_for_restart=true)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>IBFEMethod</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_f_e_method.html</anchorfile>
      <anchor>adf8197109cdbae0ac952c6b9d1aacd77</anchor>
      <arglist>(const std::string &amp;object_name, SAMRAI::tbox::Pointer&lt; SAMRAI::tbox::Database &gt; input_db, const std::vector&lt; libMesh::Mesh * &gt; &amp;meshes, int max_level_number, bool register_for_restart=true)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~IBFEMethod</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_f_e_method.html</anchorfile>
      <anchor>a6f4a24fa7e89209749ae38f386003fe1</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>IBTK::FEDataManager *</type>
      <name>getFEDataManager</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_f_e_method.html</anchorfile>
      <anchor>aa0b4c15ff591aaf6cae7426e6b8cb136</anchor>
      <arglist>(unsigned int part=0) const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>registerConstrainedPart</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_f_e_method.html</anchorfile>
      <anchor>a8304cc8f9529818f4701c831bb46c9ec</anchor>
      <arglist>(unsigned int part=0)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>registerConstrainedVelocityFunction</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_f_e_method.html</anchorfile>
      <anchor>af226af5e9fd0a23c9c621e921e5666f7</anchor>
      <arglist>(ConstrainedVelocityFcnPtr fcn, void *ctx=NULL, unsigned int part=0)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>registerConstrainedVelocityFunction</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_f_e_method.html</anchorfile>
      <anchor>aa917a579a1d4b160d68f9434db90a690</anchor>
      <arglist>(const ConstrainedVelocityFcnData &amp;data, unsigned int part=0)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>registerInitialCoordinateMappingFunction</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_f_e_method.html</anchorfile>
      <anchor>af3d9c42fed4e526885b3c763fb9cdfe4</anchor>
      <arglist>(CoordinateMappingFcnPtr fcn, void *ctx=NULL, unsigned int part=0)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>registerInitialCoordinateMappingFunction</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_f_e_method.html</anchorfile>
      <anchor>ab8c65bce8a7f44747beec40d021f0259</anchor>
      <arglist>(const CoordinateMappingFcnData &amp;data, unsigned int part=0)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>registerPK1StressFunction</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_f_e_method.html</anchorfile>
      <anchor>a88eba44b608a661eb3f3e6df255f1ef9</anchor>
      <arglist>(PK1StressFcnPtr fcn, const std::vector&lt; unsigned int &gt; &amp;systems=std::vector&lt; unsigned int &gt;(), void *ctx=NULL, libMeshEnums::QuadratureType quad_type=INVALID_Q_RULE, libMeshEnums::Order quad_order=INVALID_ORDER, unsigned int part=0)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>registerPK1StressFunction</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_f_e_method.html</anchorfile>
      <anchor>adc982545286cb33e223a44342c5455aa</anchor>
      <arglist>(const PK1StressFcnData &amp;data, unsigned int part=0)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>registerLagBodyForceFunction</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_f_e_method.html</anchorfile>
      <anchor>a36ff55504ffc10d2b409d73c15cbfc2b</anchor>
      <arglist>(LagBodyForceFcnPtr fcn, const std::vector&lt; unsigned int &gt; &amp;systems=std::vector&lt; unsigned int &gt;(), void *ctx=NULL, unsigned int part=0)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>registerLagBodyForceFunction</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_f_e_method.html</anchorfile>
      <anchor>a0e29a780e7a9973a952297317ef4537d</anchor>
      <arglist>(const LagBodyForceFcnData &amp;data, unsigned int part=0)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>registerLagSurfacePressureFunction</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_f_e_method.html</anchorfile>
      <anchor>abe317cf9e98f5754998d866c7dc773d5</anchor>
      <arglist>(LagSurfacePressureFcnPtr fcn, const std::vector&lt; unsigned int &gt; &amp;systems=std::vector&lt; unsigned int &gt;(), void *ctx=NULL, unsigned int part=0)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>registerLagSurfacePressureFunction</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_f_e_method.html</anchorfile>
      <anchor>ab62396713ffa31703aee71bc0c3e10cc</anchor>
      <arglist>(const LagSurfacePressureFcnData &amp;data, unsigned int part=0)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>registerLagSurfaceForceFunction</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_f_e_method.html</anchorfile>
      <anchor>a0bf2d6c94d6a035c6328ca4a5b716a57</anchor>
      <arglist>(LagSurfaceForceFcnPtr fcn, const std::vector&lt; unsigned int &gt; &amp;systems=std::vector&lt; unsigned int &gt;(), void *ctx=NULL, unsigned int part=0)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>registerLagSurfaceForceFunction</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_f_e_method.html</anchorfile>
      <anchor>a00593d0eebc3e8c8f4d88c99f830fd2d</anchor>
      <arglist>(const LagSurfaceForceFcnData &amp;data, unsigned int part=0)</arglist>
    </member>
    <member kind="function">
      <type>const SAMRAI::hier::IntVector&lt; NDIM &gt; &amp;</type>
      <name>getMinimumGhostCellWidth</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_f_e_method.html</anchorfile>
      <anchor>abbbb37ed429644d6f387c355d5594839</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setupTagBuffer</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_f_e_method.html</anchorfile>
      <anchor>a511797df10f977d57298d15ff1685b96</anchor>
      <arglist>(SAMRAI::tbox::Array&lt; int &gt; &amp;tag_buffer, SAMRAI::tbox::Pointer&lt; SAMRAI::mesh::GriddingAlgorithm&lt; NDIM &gt; &gt; gridding_alg) const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>preprocessIntegrateData</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_f_e_method.html</anchorfile>
      <anchor>a12bc071051c48ea46157b5cc7d7c6767</anchor>
      <arglist>(double current_time, double new_time, int num_cycles)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>postprocessIntegrateData</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_f_e_method.html</anchorfile>
      <anchor>afeb92c23596eb0fd99afe7571794642c</anchor>
      <arglist>(double current_time, double new_time, int num_cycles)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>interpolateVelocity</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_f_e_method.html</anchorfile>
      <anchor>a596f1579de307c36401a14285775809d</anchor>
      <arglist>(int u_data_idx, const std::vector&lt; SAMRAI::tbox::Pointer&lt; SAMRAI::xfer::CoarsenSchedule&lt; NDIM &gt; &gt; &gt; &amp;u_synch_scheds, const std::vector&lt; SAMRAI::tbox::Pointer&lt; SAMRAI::xfer::RefineSchedule&lt; NDIM &gt; &gt; &gt; &amp;u_ghost_fill_scheds, double data_time)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>eulerStep</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_f_e_method.html</anchorfile>
      <anchor>a8bf615ff1d61f3f769329ae1b867b19d</anchor>
      <arglist>(double current_time, double new_time)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>midpointStep</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_f_e_method.html</anchorfile>
      <anchor>a71e57708d2655dfef39026ba6bfe6b8f</anchor>
      <arglist>(double current_time, double new_time)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>trapezoidalStep</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_f_e_method.html</anchorfile>
      <anchor>a3ea509632bc801de7807aea5004d4339</anchor>
      <arglist>(double current_time, double new_time)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>computeLagrangianForce</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_f_e_method.html</anchorfile>
      <anchor>a23af664538ccd71e29ec4c03acc4a424</anchor>
      <arglist>(double data_time)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>spreadForce</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_f_e_method.html</anchorfile>
      <anchor>afa63d9f2892551528ff1eb43635e1ab1</anchor>
      <arglist>(int f_data_idx, IBTK::RobinPhysBdryPatchStrategy *f_phys_bdry_op, const std::vector&lt; SAMRAI::tbox::Pointer&lt; SAMRAI::xfer::RefineSchedule&lt; NDIM &gt; &gt; &gt; &amp;f_prolongation_scheds, double data_time)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>initializeFEData</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_f_e_method.html</anchorfile>
      <anchor>ac1fdc5d05e0f9081066b8da0be0b7e9f</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>initializePatchHierarchy</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_f_e_method.html</anchorfile>
      <anchor>a6f6d051589e682493b9fae58c129c4f3</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::hier::PatchHierarchy&lt; NDIM &gt; &gt; hierarchy, SAMRAI::tbox::Pointer&lt; SAMRAI::mesh::GriddingAlgorithm&lt; NDIM &gt; &gt; gridding_alg, int u_data_idx, const std::vector&lt; SAMRAI::tbox::Pointer&lt; SAMRAI::xfer::CoarsenSchedule&lt; NDIM &gt; &gt; &gt; &amp;u_synch_scheds, const std::vector&lt; SAMRAI::tbox::Pointer&lt; SAMRAI::xfer::RefineSchedule&lt; NDIM &gt; &gt; &gt; &amp;u_ghost_fill_scheds, int integrator_step, double init_data_time, bool initial_time)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>registerLoadBalancer</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_f_e_method.html</anchorfile>
      <anchor>ac569e03a2d23506c7563f5ec000ee594</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::mesh::LoadBalancer&lt; NDIM &gt; &gt; load_balancer, int workload_data_idx)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>updateWorkloadEstimates</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_f_e_method.html</anchorfile>
      <anchor>af7fff5d9d74aeb44d83ded26aef2a279</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::hier::PatchHierarchy&lt; NDIM &gt; &gt; hierarchy, int workload_data_idx)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>beginDataRedistribution</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_f_e_method.html</anchorfile>
      <anchor>a50a831079420eac1ccdd34c0c717933b</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::hier::PatchHierarchy&lt; NDIM &gt; &gt; hierarchy, SAMRAI::tbox::Pointer&lt; SAMRAI::mesh::GriddingAlgorithm&lt; NDIM &gt; &gt; gridding_alg)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>endDataRedistribution</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_f_e_method.html</anchorfile>
      <anchor>a8008a776c6a4f86bd516c9268c9949e4</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::hier::PatchHierarchy&lt; NDIM &gt; &gt; hierarchy, SAMRAI::tbox::Pointer&lt; SAMRAI::mesh::GriddingAlgorithm&lt; NDIM &gt; &gt; gridding_alg)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>initializeLevelData</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_f_e_method.html</anchorfile>
      <anchor>ab1c7797051fc2a1a7cac4895610d4200</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::hier::BasePatchHierarchy&lt; NDIM &gt; &gt; hierarchy, int level_number, double init_data_time, bool can_be_refined, bool initial_time, SAMRAI::tbox::Pointer&lt; SAMRAI::hier::BasePatchLevel&lt; NDIM &gt; &gt; old_level, bool allocate_data)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>resetHierarchyConfiguration</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_f_e_method.html</anchorfile>
      <anchor>aad3245b04c0f8c7b0c4066f3b099c568</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::hier::BasePatchHierarchy&lt; NDIM &gt; &gt; hierarchy, int coarsest_level, int finest_level)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>applyGradientDetector</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_f_e_method.html</anchorfile>
      <anchor>a66f3646a524f7a32cccaf62e3c370403</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::hier::BasePatchHierarchy&lt; NDIM &gt; &gt; hierarchy, int level_number, double error_data_time, int tag_index, bool initial_time, bool uses_richardson_extrapolation_too)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>putToDatabase</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_f_e_method.html</anchorfile>
      <anchor>a13fa5935694ce36da37b58ea51a4cda7</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::tbox::Database &gt; db)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>registerEulerianVariables</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_strategy.html</anchorfile>
      <anchor>ae14a796dfe21df1e9c5b5140b6573f14</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>registerEulerianCommunicationAlgorithms</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_strategy.html</anchorfile>
      <anchor>a02682292bbc009e2cd7f271a78f2c613</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>updateFixedLEOperators</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_strategy.html</anchorfile>
      <anchor>a10e59e2ec196e687c357774bb9cf2843</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual bool</type>
      <name>hasFluidSources</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_strategy.html</anchorfile>
      <anchor>aba99e5b335408d7953f4ee4340bb14e5</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>computeLagrangianFluidSource</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_strategy.html</anchorfile>
      <anchor>a9c696d55ff1c9579ac07cbca161a8b84</anchor>
      <arglist>(double data_time)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>spreadFluidSource</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_strategy.html</anchorfile>
      <anchor>a0015e5aed48df1ab5d2731e562a6b9a2</anchor>
      <arglist>(int q_data_idx, const std::vector&lt; SAMRAI::tbox::Pointer&lt; SAMRAI::xfer::RefineSchedule&lt; NDIM &gt; &gt; &gt; &amp;q_prolongation_scheds, double data_time)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>interpolatePressure</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_strategy.html</anchorfile>
      <anchor>ac2e1195ac96b299f1464d7f82f7571ac</anchor>
      <arglist>(int p_data_idx, const std::vector&lt; SAMRAI::tbox::Pointer&lt; SAMRAI::xfer::CoarsenSchedule&lt; NDIM &gt; &gt; &gt; &amp;p_synch_scheds, const std::vector&lt; SAMRAI::tbox::Pointer&lt; SAMRAI::xfer::RefineSchedule&lt; NDIM &gt; &gt; &gt; &amp;p_ghost_fill_scheds, double data_time)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>postprocessData</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_strategy.html</anchorfile>
      <anchor>a5df937f74c0871fb06fc86c284481467</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" protection="protected">
      <type>void</type>
      <name>spreadTransmissionForceDensity</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_f_e_method.html</anchorfile>
      <anchor>a67049821f3114a571a5eea4ab4c0ceaa</anchor>
      <arglist>(int f_data_idx, libMesh::PetscVector&lt; double &gt; &amp;X_ghost_vec, IBTK::RobinPhysBdryPatchStrategy *f_phys_bdry_op, double data_time, unsigned int part)</arglist>
    </member>
    <member kind="function" protection="protected">
      <type>void</type>
      <name>imposeJumpConditions</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_f_e_method.html</anchorfile>
      <anchor>af6b21ebfb951f37d4c927785db0b64ab</anchor>
      <arglist>(int f_data_idx, libMesh::PetscVector&lt; double &gt; &amp;F_ghost_vec, libMesh::PetscVector&lt; double &gt; &amp;X_ghost_vec, double data_time, unsigned int part)</arglist>
    </member>
    <member kind="function" protection="protected">
      <type>void</type>
      <name>initializeCoordinates</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_f_e_method.html</anchorfile>
      <anchor>a3f040b0cef9072dcab8874343714c851</anchor>
      <arglist>(unsigned int part)</arglist>
    </member>
    <member kind="function" protection="protected">
      <type>void</type>
      <name>updateCoordinateMapping</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_f_e_method.html</anchorfile>
      <anchor>a535ed133339bda39890af203a82d8651</anchor>
      <arglist>(unsigned int part)</arglist>
    </member>
  </compound>
  <compound kind="struct">
    <name>IBAMR::IBFEMethod::ConstrainedVelocityFcnData</name>
    <filename>struct_i_b_a_m_r_1_1_i_b_f_e_method_1_1_constrained_velocity_fcn_data.html</filename>
  </compound>
  <compound kind="struct">
    <name>IBAMR::IBFEMethod::CoordinateMappingFcnData</name>
    <filename>struct_i_b_a_m_r_1_1_i_b_f_e_method_1_1_coordinate_mapping_fcn_data.html</filename>
  </compound>
  <compound kind="struct">
    <name>IBAMR::IBFEMethod::LagBodyForceFcnData</name>
    <filename>struct_i_b_a_m_r_1_1_i_b_f_e_method_1_1_lag_body_force_fcn_data.html</filename>
  </compound>
  <compound kind="struct">
    <name>IBAMR::IBFEMethod::LagSurfaceForceFcnData</name>
    <filename>struct_i_b_a_m_r_1_1_i_b_f_e_method_1_1_lag_surface_force_fcn_data.html</filename>
  </compound>
  <compound kind="struct">
    <name>IBAMR::IBFEMethod::LagSurfacePressureFcnData</name>
    <filename>struct_i_b_a_m_r_1_1_i_b_f_e_method_1_1_lag_surface_pressure_fcn_data.html</filename>
  </compound>
  <compound kind="struct">
    <name>IBAMR::IBFEMethod::PK1StressFcnData</name>
    <filename>struct_i_b_a_m_r_1_1_i_b_f_e_method_1_1_p_k1_stress_fcn_data.html</filename>
  </compound>
  <compound kind="class">
    <name>IBAMR::IBFEPatchRecoveryPostProcessor</name>
    <filename>class_i_b_a_m_r_1_1_i_b_f_e_patch_recovery_post_processor.html</filename>
    <member kind="function">
      <type></type>
      <name>IBFEPatchRecoveryPostProcessor</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_f_e_patch_recovery_post_processor.html</anchorfile>
      <anchor>a4ee64afeedd4a9c8990447e8f9a288e4</anchor>
      <arglist>(libMesh::MeshBase *mesh, IBTK::FEDataManager *fe_data_manager)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~IBFEPatchRecoveryPostProcessor</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_f_e_patch_recovery_post_processor.html</anchorfile>
      <anchor>a1145d4946123601f310ef144516fe195</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>initializeFEData</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_f_e_patch_recovery_post_processor.html</anchorfile>
      <anchor>a60c83c118d4d5031f431aa87c56c130e</anchor>
      <arglist>(const libMesh::PeriodicBoundaries *periodic_boundaries=NULL)</arglist>
    </member>
    <member kind="function">
      <type>libMesh::System *</type>
      <name>initializeCauchyStressSystem</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_f_e_patch_recovery_post_processor.html</anchorfile>
      <anchor>ae33e084fa154d69854c4eca3e9f3c7e8</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>libMesh::System *</type>
      <name>initializePressureSystem</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_f_e_patch_recovery_post_processor.html</anchorfile>
      <anchor>a7f4e1b58f2563adc9b41931661be2bca</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>registerCauchyStressValue</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_f_e_patch_recovery_post_processor.html</anchorfile>
      <anchor>ac2c407c3ab6bddab1ec504ed7d306909</anchor>
      <arglist>(const libMesh::Elem *elem, const libMesh::QBase *qrule, unsigned int qp, const libMesh::TensorValue&lt; double &gt; &amp;sigma)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>registerPressureValue</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_f_e_patch_recovery_post_processor.html</anchorfile>
      <anchor>aa89aa93eae6bde19540887fae55cbdcd</anchor>
      <arglist>(const libMesh::Elem *elem, const libMesh::QBase *qrule, unsigned int qp, double p)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>reconstructCauchyStress</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_f_e_patch_recovery_post_processor.html</anchorfile>
      <anchor>a5f655f93fd98edca28f8e61c7ef14bef</anchor>
      <arglist>(libMesh::System &amp;sigma_system)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>reconstructPressure</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_f_e_patch_recovery_post_processor.html</anchorfile>
      <anchor>a1583a89e8180a9fbcd5347da2019a760</anchor>
      <arglist>(libMesh::System &amp;p_system)</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>IBAMR::IBFEPostProcessor</name>
    <filename>class_i_b_a_m_r_1_1_i_b_f_e_post_processor.html</filename>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>registerScalarVariable</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_f_e_post_processor.html</anchorfile>
      <anchor>af97e9a20a96425583e47589d3b361739</anchor>
      <arglist>(const std::string &amp;var_name, libMeshEnums::FEFamily var_fe_family, libMeshEnums::Order var_fe_order, IBTK::ScalarMeshFcnPtr var_fcn, std::vector&lt; unsigned int &gt; var_fcn_systems=std::vector&lt; unsigned int &gt;(), void *var_fcn_ctx=NULL)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>registerVectorVariable</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_f_e_post_processor.html</anchorfile>
      <anchor>a64f8901891943782f4a371c1bf66c776</anchor>
      <arglist>(const std::string &amp;var_name, libMeshEnums::FEFamily var_fe_family, libMeshEnums::Order var_fe_order, IBTK::VectorMeshFcnPtr var_fcn, std::vector&lt; unsigned int &gt; var_fcn_systems=std::vector&lt; unsigned int &gt;(), void *var_fcn_ctx=NULL, unsigned int var_dim=NDIM)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>registerTensorVariable</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_f_e_post_processor.html</anchorfile>
      <anchor>acdc8daf747fca4ffcece9bb54906560f</anchor>
      <arglist>(const std::string &amp;var_name, libMeshEnums::FEFamily var_fe_family, libMeshEnums::Order var_fe_order, IBTK::TensorMeshFcnPtr var_fcn, std::vector&lt; unsigned int &gt; var_fcn_systems=std::vector&lt; unsigned int &gt;(), void *var_fcn_ctx=NULL, unsigned int var_dim=NDIM)</arglist>
    </member>
    <member kind="function" protection="protected" virtualness="pure">
      <type>virtual void</type>
      <name>reconstructVariables</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_f_e_post_processor.html</anchorfile>
      <anchor>aacc05ab02693d900af693d560ae6b02c</anchor>
      <arglist>(double data_time)=0</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>IBAMR::IBHierarchyIntegrator</name>
    <filename>class_i_b_a_m_r_1_1_i_b_hierarchy_integrator.html</filename>
    <base>IBTK::HierarchyIntegrator</base>
    <class kind="class">IBAMR::IBHierarchyIntegrator::IBEulerianForceFunction</class>
    <class kind="class">IBAMR::IBHierarchyIntegrator::IBEulerianSourceFunction</class>
  </compound>
  <compound kind="class">
    <name>IBAMR::IBHierarchyIntegrator::IBEulerianForceFunction</name>
    <filename>class_i_b_a_m_r_1_1_i_b_hierarchy_integrator_1_1_i_b_eulerian_force_function.html</filename>
    <base>IBTK::CartGridFunction</base>
    <member kind="function">
      <type></type>
      <name>IBEulerianForceFunction</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_hierarchy_integrator_1_1_i_b_eulerian_force_function.html</anchorfile>
      <anchor>adcd0dbd587fe7ec45705a023a556d720</anchor>
      <arglist>(const IBHierarchyIntegrator *ib_solver)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~IBEulerianForceFunction</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_hierarchy_integrator_1_1_i_b_eulerian_force_function.html</anchorfile>
      <anchor>a1b6b05a3d6d2178e8f6faea09c382cfa</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>isTimeDependent</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_hierarchy_integrator_1_1_i_b_eulerian_force_function.html</anchorfile>
      <anchor>af4a96d8a977043ef30674cee57fa1334</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setDataOnPatchHierarchy</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_hierarchy_integrator_1_1_i_b_eulerian_force_function.html</anchorfile>
      <anchor>a84f2d4275ed3668c141a7cb525d545c3</anchor>
      <arglist>(const int data_idx, SAMRAI::tbox::Pointer&lt; SAMRAI::hier::Variable&lt; NDIM &gt; &gt; var, SAMRAI::tbox::Pointer&lt; SAMRAI::hier::PatchHierarchy&lt; NDIM &gt; &gt; hierarchy, const double data_time, const bool initial_time=false, const int coarsest_ln=-1, const int finest_ln=-1)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setDataOnPatch</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_hierarchy_integrator_1_1_i_b_eulerian_force_function.html</anchorfile>
      <anchor>af4296fc05f38befd148ef31d96ea6a84</anchor>
      <arglist>(int data_idx, SAMRAI::tbox::Pointer&lt; SAMRAI::hier::Variable&lt; NDIM &gt; &gt; var, SAMRAI::tbox::Pointer&lt; SAMRAI::hier::Patch&lt; NDIM &gt; &gt; patch, double data_time, bool initial_time=false, SAMRAI::tbox::Pointer&lt; SAMRAI::hier::PatchLevel&lt; NDIM &gt; &gt; level=SAMRAI::tbox::Pointer&lt; SAMRAI::hier::PatchLevel&lt; NDIM &gt; &gt;(NULL))</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>IBAMR::IBHierarchyIntegrator::IBEulerianSourceFunction</name>
    <filename>class_i_b_a_m_r_1_1_i_b_hierarchy_integrator_1_1_i_b_eulerian_source_function.html</filename>
    <base>IBTK::CartGridFunction</base>
    <member kind="function">
      <type></type>
      <name>IBEulerianSourceFunction</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_hierarchy_integrator_1_1_i_b_eulerian_source_function.html</anchorfile>
      <anchor>a810e57480b64777d9f570763c5788470</anchor>
      <arglist>(const IBHierarchyIntegrator *ib_solver)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~IBEulerianSourceFunction</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_hierarchy_integrator_1_1_i_b_eulerian_source_function.html</anchorfile>
      <anchor>a2207fc1d6f85fb7694b4ffb34ee0b227</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>isTimeDependent</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_hierarchy_integrator_1_1_i_b_eulerian_source_function.html</anchorfile>
      <anchor>a0cfaa3ae92622dc5fdd3766eb35bf0e3</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setDataOnPatch</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_hierarchy_integrator_1_1_i_b_eulerian_source_function.html</anchorfile>
      <anchor>a4f376303fb13cdd2bbfea2b5179a3e65</anchor>
      <arglist>(int data_idx, SAMRAI::tbox::Pointer&lt; SAMRAI::hier::Variable&lt; NDIM &gt; &gt; var, SAMRAI::tbox::Pointer&lt; SAMRAI::hier::Patch&lt; NDIM &gt; &gt; patch, double data_time, bool initial_time=false, SAMRAI::tbox::Pointer&lt; SAMRAI::hier::PatchLevel&lt; NDIM &gt; &gt; level=SAMRAI::tbox::Pointer&lt; SAMRAI::hier::PatchLevel&lt; NDIM &gt; &gt;(NULL))</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>IBAMR::IBImplicitStaggeredHierarchyIntegrator</name>
    <filename>class_i_b_a_m_r_1_1_i_b_implicit_staggered_hierarchy_integrator.html</filename>
    <base>IBAMR::IBHierarchyIntegrator</base>
    <member kind="function">
      <type></type>
      <name>IBImplicitStaggeredHierarchyIntegrator</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_implicit_staggered_hierarchy_integrator.html</anchorfile>
      <anchor>ae4ef7ee35a978db55754b4b24c86a113</anchor>
      <arglist>(const std::string &amp;object_name, SAMRAI::tbox::Pointer&lt; SAMRAI::tbox::Database &gt; input_db, SAMRAI::tbox::Pointer&lt; IBImplicitStrategy &gt; ib_method_ops, SAMRAI::tbox::Pointer&lt; INSStaggeredHierarchyIntegrator &gt; ins_hier_integrator, bool register_for_restart=true)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~IBImplicitStaggeredHierarchyIntegrator</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_implicit_staggered_hierarchy_integrator.html</anchorfile>
      <anchor>a34b2bedf163456cb2af96a8a73decfd3</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>preprocessIntegrateHierarchy</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_implicit_staggered_hierarchy_integrator.html</anchorfile>
      <anchor>a3f07f31fbc99b84cca7bcf516cac740f</anchor>
      <arglist>(double current_time, double new_time, int num_cycles=1)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>integrateHierarchy</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_implicit_staggered_hierarchy_integrator.html</anchorfile>
      <anchor>a05a288ec979541020a55a05f287f1adc</anchor>
      <arglist>(double current_time, double new_time, int cycle_num=0)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>postprocessIntegrateHierarchy</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_implicit_staggered_hierarchy_integrator.html</anchorfile>
      <anchor>a1bf1f9daf0b2d2e3f348bb0d2b92b63b</anchor>
      <arglist>(double current_time, double new_time, bool skip_synchronize_new_state_data, int num_cycles=1)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>initializeHierarchyIntegrator</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_implicit_staggered_hierarchy_integrator.html</anchorfile>
      <anchor>aed29be14f059af4cf3fdc776972d43b0</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::hier::PatchHierarchy&lt; NDIM &gt; &gt; hierarchy, SAMRAI::tbox::Pointer&lt; SAMRAI::mesh::GriddingAlgorithm&lt; NDIM &gt; &gt; gridding_alg)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>getNumberOfCycles</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_implicit_staggered_hierarchy_integrator.html</anchorfile>
      <anchor>a2aba331d9a2a0d10e41dce598c00af3f</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function" protection="protected">
      <type>void</type>
      <name>putToDatabaseSpecialized</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_implicit_staggered_hierarchy_integrator.html</anchorfile>
      <anchor>a6fcb4d9130b2d54dfc918430e6169308</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::tbox::Database &gt; db)</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>IBAMR::IBImplicitStrategy</name>
    <filename>class_i_b_a_m_r_1_1_i_b_implicit_strategy.html</filename>
    <base>IBAMR::IBStrategy</base>
    <member kind="function" virtualness="pure">
      <type>virtual void</type>
      <name>createSolverVecs</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_implicit_strategy.html</anchorfile>
      <anchor>a79240ac88fac01b7c3848b8756e205b1</anchor>
      <arglist>(Vec &amp;X_vec, Vec &amp;F_vec)=0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual void</type>
      <name>setupSolverVecs</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_implicit_strategy.html</anchorfile>
      <anchor>a1e6cc14e3d8a5fdb880fc66ca9616926</anchor>
      <arglist>(Vec &amp;X_vec, Vec &amp;F_vec)=0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual void</type>
      <name>setUpdatedPosition</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_implicit_strategy.html</anchorfile>
      <anchor>abdd827f407e056f50e665cbe59c901c8</anchor>
      <arglist>(Vec &amp;X_new_vec)=0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual void</type>
      <name>setLinearizedPosition</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_implicit_strategy.html</anchorfile>
      <anchor>a0ee7f9f5f1fe4122da9787d44e2f99b3</anchor>
      <arglist>(Vec &amp;X_vec)=0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual void</type>
      <name>computeResidual</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_implicit_strategy.html</anchorfile>
      <anchor>a61d94e4af3817e6330ece0c88b4e115e</anchor>
      <arglist>(Vec &amp;R_vec)=0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual void</type>
      <name>computeLinearizedResidual</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_implicit_strategy.html</anchorfile>
      <anchor>adda826da2e8413f53040e8cc98ddff09</anchor>
      <arglist>(Vec &amp;X_vec, Vec &amp;R_vec)=0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual void</type>
      <name>interpolateLinearizedVelocity</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_implicit_strategy.html</anchorfile>
      <anchor>a716769e239bc524b23fed5f2d0b3b0fd</anchor>
      <arglist>(int u_data_idx, const std::vector&lt; SAMRAI::tbox::Pointer&lt; SAMRAI::xfer::CoarsenSchedule&lt; NDIM &gt; &gt; &gt; &amp;u_synch_scheds, const std::vector&lt; SAMRAI::tbox::Pointer&lt; SAMRAI::xfer::RefineSchedule&lt; NDIM &gt; &gt; &gt; &amp;u_ghost_fill_scheds, double data_time)=0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual void</type>
      <name>computeLinearizedLagrangianForce</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_implicit_strategy.html</anchorfile>
      <anchor>a7ccd87e28a4347ee3ffaf2aa3771a391</anchor>
      <arglist>(Vec &amp;X_vec, double data_time)=0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual void</type>
      <name>spreadLinearizedForce</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_implicit_strategy.html</anchorfile>
      <anchor>a1cffaa80d9abce4f18b68e498266cee8</anchor>
      <arglist>(int f_data_idx, IBTK::RobinPhysBdryPatchStrategy *f_phys_bdry_op, const std::vector&lt; SAMRAI::tbox::Pointer&lt; SAMRAI::xfer::RefineSchedule&lt; NDIM &gt; &gt; &gt; &amp;f_prolongation_scheds, double data_time)=0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual const SAMRAI::hier::IntVector&lt; NDIM &gt; &amp;</type>
      <name>getMinimumGhostCellWidth</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_strategy.html</anchorfile>
      <anchor>a4eccd737f9f717fbf7b5c3c5ff87f626</anchor>
      <arglist>() const =0</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>setupTagBuffer</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_strategy.html</anchorfile>
      <anchor>a8263c7b8e9c3146aa720a1e322432b40</anchor>
      <arglist>(SAMRAI::tbox::Array&lt; int &gt; &amp;tag_buffer, SAMRAI::tbox::Pointer&lt; SAMRAI::mesh::GriddingAlgorithm&lt; NDIM &gt; &gt; gridding_alg) const </arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>preprocessIntegrateData</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_strategy.html</anchorfile>
      <anchor>a33a850e94dd2fa44854c7c99328bfb90</anchor>
      <arglist>(double current_time, double new_time, int num_cycles)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>postprocessIntegrateData</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_strategy.html</anchorfile>
      <anchor>ace0b6b640b496f34a9953569bbcb3757</anchor>
      <arglist>(double current_time, double new_time, int num_cycles)</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual void</type>
      <name>interpolateVelocity</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_strategy.html</anchorfile>
      <anchor>a5d7d2f79ab629e93a808becc43706d5b</anchor>
      <arglist>(int u_data_idx, const std::vector&lt; SAMRAI::tbox::Pointer&lt; SAMRAI::xfer::CoarsenSchedule&lt; NDIM &gt; &gt; &gt; &amp;u_synch_scheds, const std::vector&lt; SAMRAI::tbox::Pointer&lt; SAMRAI::xfer::RefineSchedule&lt; NDIM &gt; &gt; &gt; &amp;u_ghost_fill_scheds, double data_time)=0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual void</type>
      <name>eulerStep</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_strategy.html</anchorfile>
      <anchor>ab9ce696f8c76a75a44971d6bf8cf1577</anchor>
      <arglist>(double current_time, double new_time)=0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual void</type>
      <name>midpointStep</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_strategy.html</anchorfile>
      <anchor>adf7fad570f9f4757b7818ba14e5f0b66</anchor>
      <arglist>(double current_time, double new_time)=0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual void</type>
      <name>trapezoidalStep</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_strategy.html</anchorfile>
      <anchor>a426420e5ab54e4a19763eb1f13d45949</anchor>
      <arglist>(double current_time, double new_time)=0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual void</type>
      <name>computeLagrangianForce</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_strategy.html</anchorfile>
      <anchor>aec6e7fc096c898aab511e8149940618a</anchor>
      <arglist>(double data_time)=0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual void</type>
      <name>spreadForce</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_strategy.html</anchorfile>
      <anchor>a01e52f894bcd549503fe20f58dfbd08a</anchor>
      <arglist>(int f_data_idx, IBTK::RobinPhysBdryPatchStrategy *f_phys_bdry_op, const std::vector&lt; SAMRAI::tbox::Pointer&lt; SAMRAI::xfer::RefineSchedule&lt; NDIM &gt; &gt; &gt; &amp;f_prolongation_scheds, double data_time)=0</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>initializePatchHierarchy</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_strategy.html</anchorfile>
      <anchor>a7fbf36ad9f32795bd640ce2bab81f93a</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::hier::PatchHierarchy&lt; NDIM &gt; &gt; hierarchy, SAMRAI::tbox::Pointer&lt; SAMRAI::mesh::GriddingAlgorithm&lt; NDIM &gt; &gt; gridding_alg, int u_data_idx, const std::vector&lt; SAMRAI::tbox::Pointer&lt; SAMRAI::xfer::CoarsenSchedule&lt; NDIM &gt; &gt; &gt; &amp;u_synch_scheds, const std::vector&lt; SAMRAI::tbox::Pointer&lt; SAMRAI::xfer::RefineSchedule&lt; NDIM &gt; &gt; &gt; &amp;u_ghost_fill_scheds, int integrator_step, double init_data_time, bool initial_time)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>registerLoadBalancer</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_strategy.html</anchorfile>
      <anchor>a2d7cb774d87cdd04e2ccda9f4f7daddb</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::mesh::LoadBalancer&lt; NDIM &gt; &gt; load_balancer, int workload_data_idx)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>updateWorkloadEstimates</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_strategy.html</anchorfile>
      <anchor>a5e2fb29d1a2278cf8dadf7da019a26c2</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::hier::PatchHierarchy&lt; NDIM &gt; &gt; hierarchy, int workload_data_idx)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>beginDataRedistribution</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_strategy.html</anchorfile>
      <anchor>a831d4e247e9b641d6b11f67b686924de</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::hier::PatchHierarchy&lt; NDIM &gt; &gt; hierarchy, SAMRAI::tbox::Pointer&lt; SAMRAI::mesh::GriddingAlgorithm&lt; NDIM &gt; &gt; gridding_alg)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>endDataRedistribution</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_strategy.html</anchorfile>
      <anchor>ab32e858a47b116152cf87a7602488796</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::hier::PatchHierarchy&lt; NDIM &gt; &gt; hierarchy, SAMRAI::tbox::Pointer&lt; SAMRAI::mesh::GriddingAlgorithm&lt; NDIM &gt; &gt; gridding_alg)</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>IBAMR::IBInstrumentationSpec</name>
    <filename>class_i_b_a_m_r_1_1_i_b_instrumentation_spec.html</filename>
    <base>IBTK::Streamable</base>
    <member kind="function">
      <type></type>
      <name>IBInstrumentationSpec</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_instrumentation_spec.html</anchorfile>
      <anchor>a68611e3c3f2f4cf1143f036344ed644d</anchor>
      <arglist>(int master_idx=-1, int meter_idx=-1, int node_idx=-1)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~IBInstrumentationSpec</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_instrumentation_spec.html</anchorfile>
      <anchor>a535e516b8c6bd4d0d225ad95753df8d9</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>const int &amp;</type>
      <name>getMasterNodeIndex</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_instrumentation_spec.html</anchorfile>
      <anchor>ac76b354f3949af7fa77a79860d112ed6</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>int &amp;</type>
      <name>getMasterNodeIndex</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_instrumentation_spec.html</anchorfile>
      <anchor>ad9ec55bfe5f4d694f356d2d58bb91d61</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>const int &amp;</type>
      <name>getMeterIndex</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_instrumentation_spec.html</anchorfile>
      <anchor>a966131547b133b018f2f84b4b85cc32b</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>int &amp;</type>
      <name>getMeterIndex</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_instrumentation_spec.html</anchorfile>
      <anchor>aa46b56292bd09ee8d895f22a01dd1202</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>const int &amp;</type>
      <name>getNodeIndex</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_instrumentation_spec.html</anchorfile>
      <anchor>af9933c8158c064c5d1636ee9af131285</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>int &amp;</type>
      <name>getNodeIndex</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_instrumentation_spec.html</anchorfile>
      <anchor>af04ad5a3ab881a8d25d4e2b621e26b62</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>getStreamableClassID</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_instrumentation_spec.html</anchorfile>
      <anchor>a6d7b59612eced7f29cd128461aa1544f</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>size_t</type>
      <name>getDataStreamSize</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_instrumentation_spec.html</anchorfile>
      <anchor>ac8a1ce182b1ac56d190b541607b0876e</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>packStream</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_instrumentation_spec.html</anchorfile>
      <anchor>aef577d8d9fbd98dada1ca94ede16e75d</anchor>
      <arglist>(SAMRAI::tbox::AbstractStream &amp;stream)</arglist>
    </member>
    <member kind="function" static="yes">
      <type>static void</type>
      <name>registerWithStreamableManager</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_instrumentation_spec.html</anchorfile>
      <anchor>a1d8b9983d77460e7971206eafd537c91</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" static="yes">
      <type>static bool</type>
      <name>getIsRegisteredWithStreamableManager</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_instrumentation_spec.html</anchorfile>
      <anchor>a09cd34b435143d3c33e62217c7df7b5a</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" static="yes">
      <type>static void</type>
      <name>setInstrumentNames</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_instrumentation_spec.html</anchorfile>
      <anchor>a376ef52a32a8116130fcefc73a703ba5</anchor>
      <arglist>(const std::vector&lt; std::string &gt; &amp;names)</arglist>
    </member>
    <member kind="function" static="yes">
      <type>static const std::vector&lt; std::string &gt; &amp;</type>
      <name>getInstrumentNames</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_instrumentation_spec.html</anchorfile>
      <anchor>a9e5ce344359123574ba1664c2ae67656</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="variable" static="yes">
      <type>static int</type>
      <name>STREAMABLE_CLASS_ID</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_instrumentation_spec.html</anchorfile>
      <anchor>aebf96a76b2d19cb1dd70575e4735785a</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>IBAMR::IBInstrumentPanel</name>
    <filename>class_i_b_a_m_r_1_1_i_b_instrument_panel.html</filename>
    <member kind="function">
      <type></type>
      <name>IBInstrumentPanel</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_instrument_panel.html</anchorfile>
      <anchor>a5df6bb5ab0b5b8472d0a18e404f0988d</anchor>
      <arglist>(const std::string &amp;object_name, SAMRAI::tbox::Pointer&lt; SAMRAI::tbox::Database &gt; input_db)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~IBInstrumentPanel</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_instrument_panel.html</anchorfile>
      <anchor>a3a2afe4f2ea6a804bdf970b47cba68b7</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>const std::vector&lt; std::string &gt; &amp;</type>
      <name>getInstrumentNames</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_instrument_panel.html</anchorfile>
      <anchor>ad71246f4e609fcbfc5940ea02914c757</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>const double &amp;</type>
      <name>getInstrumentDataReadTime</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_instrument_panel.html</anchorfile>
      <anchor>ac13bc69c5c7f81cb8dae90bce089d3b2</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>const std::vector&lt; double &gt; &amp;</type>
      <name>getFlowValues</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_instrument_panel.html</anchorfile>
      <anchor>ac8a94466f0924a3b4dc82b606413e8c9</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>const std::vector&lt; double &gt; &amp;</type>
      <name>getMeanPressureValues</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_instrument_panel.html</anchorfile>
      <anchor>a8ae77c682dc3020330a3ed2a301b1264</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>const std::vector&lt; double &gt; &amp;</type>
      <name>getPointwisePressureValues</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_instrument_panel.html</anchorfile>
      <anchor>af4d25de39a0c176b86cef25a3eb3068a</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>isInstrumented</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_instrument_panel.html</anchorfile>
      <anchor>ac167a7dc60efce26597fca19b57a0351</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>initializeHierarchyIndependentData</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_instrument_panel.html</anchorfile>
      <anchor>a830cb7ab85ab3baf80c2e31ca835e1d1</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::hier::PatchHierarchy&lt; NDIM &gt; &gt; hierarchy, IBTK::LDataManager *l_data_manager)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>initializeHierarchyDependentData</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_instrument_panel.html</anchorfile>
      <anchor>a84f8fc870fb48fdf1f826fdfc8ac6ce1</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::hier::PatchHierarchy&lt; NDIM &gt; &gt; hierarchy, IBTK::LDataManager *l_data_manager, int timestep_num, double data_time)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>readInstrumentData</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_instrument_panel.html</anchorfile>
      <anchor>a6c2d52712090e3334afb4e26f31b8e64</anchor>
      <arglist>(int U_data_idx, int P_data_idx, SAMRAI::tbox::Pointer&lt; SAMRAI::hier::PatchHierarchy&lt; NDIM &gt; &gt; hierarchy, IBTK::LDataManager *l_data_manager, int timestep_num, double data_time)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setPlotDirectory</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_instrument_panel.html</anchorfile>
      <anchor>ae2fe5b17ffda5c01e0112bbd76d89d8f</anchor>
      <arglist>(const std::string &amp;plot_directory_name)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>writePlotData</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_instrument_panel.html</anchorfile>
      <anchor>aabfe0714bac4e83ff2200ce1ef9ad05f</anchor>
      <arglist>(int timestep_num, double simulation_time)</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>IBAMR::IBKirchhoffRodForceGen</name>
    <filename>class_i_b_a_m_r_1_1_i_b_kirchhoff_rod_force_gen.html</filename>
    <member kind="function">
      <type></type>
      <name>IBKirchhoffRodForceGen</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_kirchhoff_rod_force_gen.html</anchorfile>
      <anchor>a8f56a98eca9e26713e97648c521bb83f</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::tbox::Database &gt; input_db=NULL)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~IBKirchhoffRodForceGen</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_kirchhoff_rod_force_gen.html</anchorfile>
      <anchor>ae24c580ab3496aa200b6af694c95e0f5</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>initializeLevelData</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_kirchhoff_rod_force_gen.html</anchorfile>
      <anchor>ae1f6504913882f2fad815e95a80a88ca</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::hier::PatchHierarchy&lt; NDIM &gt; &gt; hierarchy, int level_number, double init_data_time, bool initial_time, IBTK::LDataManager *l_data_manager)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>computeLagrangianForceAndTorque</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_kirchhoff_rod_force_gen.html</anchorfile>
      <anchor>a12e027be6800df9bf373f829fac99efe</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; IBTK::LData &gt; F_data, SAMRAI::tbox::Pointer&lt; IBTK::LData &gt; N_data, SAMRAI::tbox::Pointer&lt; IBTK::LData &gt; X_data, SAMRAI::tbox::Pointer&lt; IBTK::LData &gt; D_data, SAMRAI::tbox::Pointer&lt; SAMRAI::hier::PatchHierarchy&lt; NDIM &gt; &gt; hierarchy, int level_number, double data_time, IBTK::LDataManager *l_data_manager)</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>IBAMR::IBLagrangianForceStrategy</name>
    <filename>class_i_b_a_m_r_1_1_i_b_lagrangian_force_strategy.html</filename>
    <member kind="function">
      <type></type>
      <name>IBLagrangianForceStrategy</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_lagrangian_force_strategy.html</anchorfile>
      <anchor>adddcaca443f205a3efef3fa2f5ae139c</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual</type>
      <name>~IBLagrangianForceStrategy</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_lagrangian_force_strategy.html</anchorfile>
      <anchor>a189300a0b465d6f4c053f3d28851be75</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>setTimeInterval</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_lagrangian_force_strategy.html</anchorfile>
      <anchor>a61bd5982b47d434a5854a0ed552ca62d</anchor>
      <arglist>(double current_time, double new_time)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>initializeLevelData</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_lagrangian_force_strategy.html</anchorfile>
      <anchor>acbd1b66f4a92a87cb349f837e327f9d0</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::hier::PatchHierarchy&lt; NDIM &gt; &gt; hierarchy, int level_number, double init_data_time, bool initial_time, IBTK::LDataManager *l_data_manager)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>computeLagrangianForce</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_lagrangian_force_strategy.html</anchorfile>
      <anchor>aadd823e0da9a01a6d24facb1db7e1922</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; IBTK::LData &gt; F_data, SAMRAI::tbox::Pointer&lt; IBTK::LData &gt; X_data, SAMRAI::tbox::Pointer&lt; IBTK::LData &gt; U_data, SAMRAI::tbox::Pointer&lt; SAMRAI::hier::PatchHierarchy&lt; NDIM &gt; &gt; hierarchy, int level_number, double data_time, IBTK::LDataManager *l_data_manager)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>computeLagrangianForceJacobianNonzeroStructure</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_lagrangian_force_strategy.html</anchorfile>
      <anchor>a8ca2cdee686ac794f83095ae17947d1f</anchor>
      <arglist>(std::vector&lt; int &gt; &amp;d_nnz, std::vector&lt; int &gt; &amp;o_nnz, SAMRAI::tbox::Pointer&lt; SAMRAI::hier::PatchHierarchy&lt; NDIM &gt; &gt; hierarchy, int level_number, IBTK::LDataManager *l_data_manager)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>computeLagrangianForceJacobian</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_lagrangian_force_strategy.html</anchorfile>
      <anchor>af2ae6cb93a2e2038060a507ab2cfa6c0</anchor>
      <arglist>(Mat &amp;J_mat, MatAssemblyType assembly_type, double X_coef, SAMRAI::tbox::Pointer&lt; IBTK::LData &gt; X_data, double U_coef, SAMRAI::tbox::Pointer&lt; IBTK::LData &gt; U_data, SAMRAI::tbox::Pointer&lt; SAMRAI::hier::PatchHierarchy&lt; NDIM &gt; &gt; hierarchy, int level_number, double data_time, IBTK::LDataManager *l_data_manager)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual double</type>
      <name>computeLagrangianEnergy</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_lagrangian_force_strategy.html</anchorfile>
      <anchor>afdccf23d3d07ff142450f8703b79b6bf</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; IBTK::LData &gt; X_data, SAMRAI::tbox::Pointer&lt; IBTK::LData &gt; U_data, SAMRAI::tbox::Pointer&lt; SAMRAI::hier::PatchHierarchy&lt; NDIM &gt; &gt; hierarchy, int level_number, double data_time, IBTK::LDataManager *l_data_manager)</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>IBAMR::IBLagrangianForceStrategySet</name>
    <filename>class_i_b_a_m_r_1_1_i_b_lagrangian_force_strategy_set.html</filename>
    <base>IBAMR::IBLagrangianForceStrategy</base>
    <member kind="function">
      <type></type>
      <name>IBLagrangianForceStrategySet</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_lagrangian_force_strategy_set.html</anchorfile>
      <anchor>a9c144b8acd6b1e0a2aecd6a0a599e287</anchor>
      <arglist>(InputIterator first, InputIterator last)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~IBLagrangianForceStrategySet</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_lagrangian_force_strategy_set.html</anchorfile>
      <anchor>a415e45b4916fd188614eb67985d69706</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setTimeInterval</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_lagrangian_force_strategy_set.html</anchorfile>
      <anchor>a52f9b53b015b14ff63f82798679efe49</anchor>
      <arglist>(double current_time, double new_time)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>initializeLevelData</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_lagrangian_force_strategy_set.html</anchorfile>
      <anchor>ad316c06140050092208a13fe6a0bf67c</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::hier::PatchHierarchy&lt; NDIM &gt; &gt; hierarchy, int level_number, double init_data_time, bool initial_time, IBTK::LDataManager *l_data_manager)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>computeLagrangianForce</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_lagrangian_force_strategy_set.html</anchorfile>
      <anchor>a63629caed18b778ffcb4a2e2678bcebc</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; IBTK::LData &gt; F_data, SAMRAI::tbox::Pointer&lt; IBTK::LData &gt; X_data, SAMRAI::tbox::Pointer&lt; IBTK::LData &gt; U_data, SAMRAI::tbox::Pointer&lt; SAMRAI::hier::PatchHierarchy&lt; NDIM &gt; &gt; hierarchy, int level_number, double data_time, IBTK::LDataManager *l_data_manager)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>computeLagrangianForceJacobianNonzeroStructure</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_lagrangian_force_strategy_set.html</anchorfile>
      <anchor>ab2ab1729d20218f93e18aa9d0bbb208f</anchor>
      <arglist>(std::vector&lt; int &gt; &amp;d_nnz, std::vector&lt; int &gt; &amp;o_nnz, SAMRAI::tbox::Pointer&lt; SAMRAI::hier::PatchHierarchy&lt; NDIM &gt; &gt; hierarchy, int level_number, IBTK::LDataManager *l_data_manager)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>computeLagrangianForceJacobian</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_lagrangian_force_strategy_set.html</anchorfile>
      <anchor>a624e5430a3a96e156a46a0583c4a3c9a</anchor>
      <arglist>(Mat &amp;J_mat, MatAssemblyType assembly_type, double X_coef, SAMRAI::tbox::Pointer&lt; IBTK::LData &gt; X_data, double U_coef, SAMRAI::tbox::Pointer&lt; IBTK::LData &gt; U_data, SAMRAI::tbox::Pointer&lt; SAMRAI::hier::PatchHierarchy&lt; NDIM &gt; &gt; hierarchy, int level_number, double data_time, IBTK::LDataManager *l_data_manager)</arglist>
    </member>
    <member kind="function">
      <type>double</type>
      <name>computeLagrangianEnergy</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_lagrangian_force_strategy_set.html</anchorfile>
      <anchor>a10f03994ed984851b4a36d2ef3861e17</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; IBTK::LData &gt; X_data, SAMRAI::tbox::Pointer&lt; IBTK::LData &gt; U_data, SAMRAI::tbox::Pointer&lt; SAMRAI::hier::PatchHierarchy&lt; NDIM &gt; &gt; hierarchy, int level_number, double data_time, IBTK::LDataManager *l_data_manager)</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>IBAMR::IBLagrangianSourceStrategy</name>
    <filename>class_i_b_a_m_r_1_1_i_b_lagrangian_source_strategy.html</filename>
    <member kind="function">
      <type></type>
      <name>IBLagrangianSourceStrategy</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_lagrangian_source_strategy.html</anchorfile>
      <anchor>a9c1bb0e89ee7deed2c2dd5d077b767c6</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual</type>
      <name>~IBLagrangianSourceStrategy</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_lagrangian_source_strategy.html</anchorfile>
      <anchor>a4397651de8840a08ab3bd5354d474cf5</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>setTimeInterval</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_lagrangian_source_strategy.html</anchorfile>
      <anchor>af0bfa2f590cfeef97b561c9e73bbf3a3</anchor>
      <arglist>(double current_time, double new_time)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>initializeLevelData</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_lagrangian_source_strategy.html</anchorfile>
      <anchor>ab4508d9f4d72a1f2c1f6b94bf9814758</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::hier::PatchHierarchy&lt; NDIM &gt; &gt; hierarchy, int level_number, double init_data_time, bool initial_time, IBTK::LDataManager *l_data_manager)</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual unsigned int</type>
      <name>getNumSources</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_lagrangian_source_strategy.html</anchorfile>
      <anchor>a08fedc07596a4d61263dc435f4b5c7ed</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::hier::PatchHierarchy&lt; NDIM &gt; &gt; hierarchy, int level_number, double data_time, IBTK::LDataManager *l_data_manager)=0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual void</type>
      <name>getSourceLocations</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_lagrangian_source_strategy.html</anchorfile>
      <anchor>a2f0bbeeff91ed5e18ae8dbeabb51dde6</anchor>
      <arglist>(std::vector&lt; IBTK::Point &gt; &amp;X_src, std::vector&lt; double &gt; &amp;r_src, SAMRAI::tbox::Pointer&lt; IBTK::LData &gt; X_data, SAMRAI::tbox::Pointer&lt; SAMRAI::hier::PatchHierarchy&lt; NDIM &gt; &gt; hierarchy, int level_number, double data_time, IBTK::LDataManager *l_data_manager)=0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual void</type>
      <name>setSourcePressures</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_lagrangian_source_strategy.html</anchorfile>
      <anchor>aba52143da123a93f46b9edf1f3ea1ae0</anchor>
      <arglist>(const std::vector&lt; double &gt; &amp;P_src, SAMRAI::tbox::Pointer&lt; SAMRAI::hier::PatchHierarchy&lt; NDIM &gt; &gt; hierarchy, int level_number, double data_time, IBTK::LDataManager *l_data_manager)=0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual void</type>
      <name>computeSourceStrengths</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_lagrangian_source_strategy.html</anchorfile>
      <anchor>ac29365ac89304816f0c2ed85cafc30e5</anchor>
      <arglist>(std::vector&lt; double &gt; &amp;Q_src, SAMRAI::tbox::Pointer&lt; SAMRAI::hier::PatchHierarchy&lt; NDIM &gt; &gt; hierarchy, int level_number, double data_time, IBTK::LDataManager *l_data_manager)=0</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>IBAMR::IBMethod</name>
    <filename>class_i_b_a_m_r_1_1_i_b_method.html</filename>
    <base>IBAMR::IBImplicitStrategy</base>
  </compound>
  <compound kind="class">
    <name>IBAMR::IBMethodPostProcessStrategy</name>
    <filename>class_i_b_a_m_r_1_1_i_b_method_post_process_strategy.html</filename>
    <member kind="function">
      <type></type>
      <name>IBMethodPostProcessStrategy</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_method_post_process_strategy.html</anchorfile>
      <anchor>aa902c325da8347b9ced62fd15e540c46</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual</type>
      <name>~IBMethodPostProcessStrategy</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_method_post_process_strategy.html</anchorfile>
      <anchor>aced5db0524a8dc26754e23ec6d25ba25</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual void</type>
      <name>postprocessData</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_method_post_process_strategy.html</anchorfile>
      <anchor>ab7ea2df0b36d87d28ce1af0d3496e560</anchor>
      <arglist>(int u_idx, int p_idx, int f_idx, const std::vector&lt; SAMRAI::tbox::Pointer&lt; IBTK::LData &gt; &gt; &amp;F_data, const std::vector&lt; SAMRAI::tbox::Pointer&lt; IBTK::LData &gt; &gt; &amp;X_data, const std::vector&lt; SAMRAI::tbox::Pointer&lt; IBTK::LData &gt; &gt; &amp;U_data, SAMRAI::tbox::Pointer&lt; SAMRAI::hier::PatchHierarchy&lt; NDIM &gt; &gt; hierarchy, int coarsest_level_number, int finest_level_number, double data_time, IBMethod *ib_method)=0</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>IBAMR::IBRodForceSpec</name>
    <filename>class_i_b_a_m_r_1_1_i_b_rod_force_spec.html</filename>
    <base>IBTK::Streamable</base>
    <member kind="function">
      <type></type>
      <name>IBRodForceSpec</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_rod_force_spec.html</anchorfile>
      <anchor>a0b4c1dc5a35020f9c76305e14326fb3a</anchor>
      <arglist>(unsigned int num_rods=0)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>IBRodForceSpec</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_rod_force_spec.html</anchorfile>
      <anchor>abaa092081e221ebee5267176da4d9c97</anchor>
      <arglist>(int master_idx, const std::vector&lt; int &gt; &amp;next_idxs, const std::vector&lt; boost::array&lt; double, NUM_MATERIAL_PARAMS &gt; &gt; &amp;material_params)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~IBRodForceSpec</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_rod_force_spec.html</anchorfile>
      <anchor>a5e8a322f7310fb1c9b03ff272d82e487</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>unsigned int</type>
      <name>getNumberOfRods</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_rod_force_spec.html</anchorfile>
      <anchor>a0f60b0f5228607a4969dfd749e168cbc</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>const int &amp;</type>
      <name>getMasterNodeIndex</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_rod_force_spec.html</anchorfile>
      <anchor>a85190c8db2ec9195d323fa402a35f353</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>int &amp;</type>
      <name>getMasterNodeIndex</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_rod_force_spec.html</anchorfile>
      <anchor>a2f9655a4e544d23b91a3d1c889566b5d</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>const std::vector&lt; int &gt; &amp;</type>
      <name>getNextNodeIndices</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_rod_force_spec.html</anchorfile>
      <anchor>aeac84e2121b20e92f8e5d47141dab19f</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>std::vector&lt; int &gt; &amp;</type>
      <name>getNextNodeIndices</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_rod_force_spec.html</anchorfile>
      <anchor>af63276865a39b3e5b0d2db42d7144da9</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>const std::vector&lt; boost::array&lt; double, NUM_MATERIAL_PARAMS &gt; &gt; &amp;</type>
      <name>getMaterialParams</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_rod_force_spec.html</anchorfile>
      <anchor>a6196ce11ec4078e54c1683cc57961af6</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>std::vector&lt; boost::array&lt; double, NUM_MATERIAL_PARAMS &gt; &gt; &amp;</type>
      <name>getMaterialParams</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_rod_force_spec.html</anchorfile>
      <anchor>a5bcf6db7307fd575fca9571850937276</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>getStreamableClassID</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_rod_force_spec.html</anchorfile>
      <anchor>aa3e2ae9c2f61377e8102e9e029e0bc8d</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>size_t</type>
      <name>getDataStreamSize</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_rod_force_spec.html</anchorfile>
      <anchor>a238a6d1557c158497cc46d46ec79e922</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>packStream</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_rod_force_spec.html</anchorfile>
      <anchor>a2fbf03ec417a83ee8fd18f1b451bbcf6</anchor>
      <arglist>(SAMRAI::tbox::AbstractStream &amp;stream)</arglist>
    </member>
    <member kind="function" static="yes">
      <type>static void</type>
      <name>registerWithStreamableManager</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_rod_force_spec.html</anchorfile>
      <anchor>acc98e9054aff777abcebd940e9d15195</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" static="yes">
      <type>static bool</type>
      <name>getIsRegisteredWithStreamableManager</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_rod_force_spec.html</anchorfile>
      <anchor>af3f0226f2a788c84c5eb469fdf5f773c</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="variable" static="yes">
      <type>static int</type>
      <name>STREAMABLE_CLASS_ID</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_rod_force_spec.html</anchorfile>
      <anchor>aa268cc812174fdc31ca2b37577f5fc5d</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>IBAMR::IBSourceSpec</name>
    <filename>class_i_b_a_m_r_1_1_i_b_source_spec.html</filename>
    <base>IBTK::Streamable</base>
    <member kind="function">
      <type></type>
      <name>IBSourceSpec</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_source_spec.html</anchorfile>
      <anchor>a0246347a6890a65fc77ae1141f3a9f0d</anchor>
      <arglist>(int master_idx=-1, int source_idx=-1)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~IBSourceSpec</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_source_spec.html</anchorfile>
      <anchor>a71a996a7ef9758e08aae36e57e6e17f4</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>const int &amp;</type>
      <name>getMasterNodeIndex</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_source_spec.html</anchorfile>
      <anchor>a7aca9b5013183144fee27939c99a3055</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>int &amp;</type>
      <name>getMasterNodeIndex</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_source_spec.html</anchorfile>
      <anchor>a5c08f0ca9f6c50b5f6bc15d43ca236f9</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>const int &amp;</type>
      <name>getSourceIndex</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_source_spec.html</anchorfile>
      <anchor>a9fad0f302ac123972731d23b00779d42</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>int &amp;</type>
      <name>getSourceIndex</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_source_spec.html</anchorfile>
      <anchor>a62534e240c13b628d914df9715fae5b9</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>getStreamableClassID</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_source_spec.html</anchorfile>
      <anchor>a274b469597a7e457e837b7207a524de9</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>size_t</type>
      <name>getDataStreamSize</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_source_spec.html</anchorfile>
      <anchor>a63666c005e53bb2a4ec18b6af618bfe1</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>packStream</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_source_spec.html</anchorfile>
      <anchor>ab980cae2507aee34e7c2bf3156f50be7</anchor>
      <arglist>(SAMRAI::tbox::AbstractStream &amp;stream)</arglist>
    </member>
    <member kind="function" static="yes">
      <type>static void</type>
      <name>registerWithStreamableManager</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_source_spec.html</anchorfile>
      <anchor>a48b32905ba75770566027a62c8e70bb6</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" static="yes">
      <type>static bool</type>
      <name>getIsRegisteredWithStreamableManager</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_source_spec.html</anchorfile>
      <anchor>a385bd7c52925887a5bbf24a1298c6006</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="variable" static="yes">
      <type>static int</type>
      <name>STREAMABLE_CLASS_ID</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_source_spec.html</anchorfile>
      <anchor>aad56f404f0164817f1cf7780746a4bee</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>IBAMR::IBSpringForceSpec</name>
    <filename>class_i_b_a_m_r_1_1_i_b_spring_force_spec.html</filename>
    <base>IBTK::Streamable</base>
    <member kind="function">
      <type></type>
      <name>IBSpringForceSpec</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_spring_force_spec.html</anchorfile>
      <anchor>a0023bb0b66b793f40a5a5568010c864d</anchor>
      <arglist>(unsigned int num_springs=0)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>IBSpringForceSpec</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_spring_force_spec.html</anchorfile>
      <anchor>a28f0c799e8348d44d52cd2242f2a302c</anchor>
      <arglist>(int master_idx, const std::vector&lt; int &gt; &amp;slave_idxs, const std::vector&lt; int &gt; &amp;force_fcn_idxs, const std::vector&lt; std::vector&lt; double &gt; &gt; &amp;parameters)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~IBSpringForceSpec</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_spring_force_spec.html</anchorfile>
      <anchor>a76f0987c3ea21fc8c505c3cb73517cfb</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>unsigned int</type>
      <name>getNumberOfSprings</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_spring_force_spec.html</anchorfile>
      <anchor>afe650d855fd2c25305ea3762d3c2e89b</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>const int &amp;</type>
      <name>getMasterNodeIndex</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_spring_force_spec.html</anchorfile>
      <anchor>afb72be1910296af50a48115f49eba903</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>int &amp;</type>
      <name>getMasterNodeIndex</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_spring_force_spec.html</anchorfile>
      <anchor>abb241b8048707ad1310f5e56d53500b2</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>const std::vector&lt; int &gt; &amp;</type>
      <name>getSlaveNodeIndices</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_spring_force_spec.html</anchorfile>
      <anchor>a607c64d4e6feaba9135816848d80c7fa</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>std::vector&lt; int &gt; &amp;</type>
      <name>getSlaveNodeIndices</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_spring_force_spec.html</anchorfile>
      <anchor>a4370c292c88fc1eb851f262a71ffaef8</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>const std::vector&lt; int &gt; &amp;</type>
      <name>getForceFunctionIndices</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_spring_force_spec.html</anchorfile>
      <anchor>ab89cd7ba2abb90d969699ecc9f55a7ed</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>std::vector&lt; int &gt; &amp;</type>
      <name>getForceFunctionIndices</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_spring_force_spec.html</anchorfile>
      <anchor>a2b4d1e7f676a24728dc6f7eec3e9ef18</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>const std::vector&lt; std::vector&lt; double &gt; &gt; &amp;</type>
      <name>getParameters</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_spring_force_spec.html</anchorfile>
      <anchor>a2050a18239bac574f3bccbb6f7cfd0d6</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>std::vector&lt; std::vector&lt; double &gt; &gt; &amp;</type>
      <name>getParameters</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_spring_force_spec.html</anchorfile>
      <anchor>adc0d8bacc65981a66729455d82d60eb0</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>getStreamableClassID</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_spring_force_spec.html</anchorfile>
      <anchor>a54d1739a8ca01f2f245a74f9c95b3583</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>size_t</type>
      <name>getDataStreamSize</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_spring_force_spec.html</anchorfile>
      <anchor>a5022d47c827b29e8aca2ed5716d99dc3</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>packStream</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_spring_force_spec.html</anchorfile>
      <anchor>aa8b966ee7511eab3527da79651e99ad4</anchor>
      <arglist>(SAMRAI::tbox::AbstractStream &amp;stream)</arglist>
    </member>
    <member kind="function" static="yes">
      <type>static void</type>
      <name>registerWithStreamableManager</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_spring_force_spec.html</anchorfile>
      <anchor>ac3bf14acfd7723b3fb42b420a6f74c50</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" static="yes">
      <type>static bool</type>
      <name>getIsRegisteredWithStreamableManager</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_spring_force_spec.html</anchorfile>
      <anchor>ac45711834fe2480ac48cd08399e327b8</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="variable" static="yes">
      <type>static int</type>
      <name>STREAMABLE_CLASS_ID</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_spring_force_spec.html</anchorfile>
      <anchor>ac0dd2ecec20304b24fe3190f9a9a90be</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>IBAMR::IBStandardForceGen</name>
    <filename>class_i_b_a_m_r_1_1_i_b_standard_force_gen.html</filename>
    <base>IBAMR::IBLagrangianForceStrategy</base>
    <member kind="function">
      <type></type>
      <name>IBStandardForceGen</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_standard_force_gen.html</anchorfile>
      <anchor>a0bf0993a80ce2519f1f454b97e10529c</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~IBStandardForceGen</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_standard_force_gen.html</anchorfile>
      <anchor>a17c73a33ad726756a984102da10ca272</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>registerSpringForceFunction</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_standard_force_gen.html</anchorfile>
      <anchor>a4ea87a5ae6b435f0503cbd0825b28530</anchor>
      <arglist>(int force_fcn_index, const SpringForceFcnPtr spring_force_fcn_ptr, const SpringForceDerivFcnPtr spring_force_deriv_fcn_ptr=NULL)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>initializeLevelData</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_standard_force_gen.html</anchorfile>
      <anchor>a46478b2413e332890e0d87f7f309a779</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::hier::PatchHierarchy&lt; NDIM &gt; &gt; hierarchy, int level_number, double init_data_time, bool initial_time, IBTK::LDataManager *l_data_manager)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>computeLagrangianForce</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_standard_force_gen.html</anchorfile>
      <anchor>af5adb349135b47911d04a85a6511b493</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; IBTK::LData &gt; F_data, SAMRAI::tbox::Pointer&lt; IBTK::LData &gt; X_data, SAMRAI::tbox::Pointer&lt; IBTK::LData &gt; U_data, SAMRAI::tbox::Pointer&lt; SAMRAI::hier::PatchHierarchy&lt; NDIM &gt; &gt; hierarchy, int level_number, double data_time, IBTK::LDataManager *l_data_manager)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>computeLagrangianForceJacobianNonzeroStructure</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_standard_force_gen.html</anchorfile>
      <anchor>aa6a48b63acd8fac23364b182a4c0b776</anchor>
      <arglist>(std::vector&lt; int &gt; &amp;d_nnz, std::vector&lt; int &gt; &amp;o_nnz, SAMRAI::tbox::Pointer&lt; SAMRAI::hier::PatchHierarchy&lt; NDIM &gt; &gt; hierarchy, int level_number, IBTK::LDataManager *l_data_manager)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>computeLagrangianForceJacobian</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_standard_force_gen.html</anchorfile>
      <anchor>aa98b660cb901f23fb8125e317c06aacf</anchor>
      <arglist>(Mat &amp;J_mat, MatAssemblyType assembly_type, double X_coef, SAMRAI::tbox::Pointer&lt; IBTK::LData &gt; X_data, double U_coef, SAMRAI::tbox::Pointer&lt; IBTK::LData &gt; U_data, SAMRAI::tbox::Pointer&lt; SAMRAI::hier::PatchHierarchy&lt; NDIM &gt; &gt; hierarchy, int level_number, double data_time, IBTK::LDataManager *l_data_manager)</arglist>
    </member>
    <member kind="function">
      <type>double</type>
      <name>computeLagrangianEnergy</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_standard_force_gen.html</anchorfile>
      <anchor>a0269c77dad917a2c18b7dbd3fa4fcbb8</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; IBTK::LData &gt; X_data, SAMRAI::tbox::Pointer&lt; IBTK::LData &gt; U_data, SAMRAI::tbox::Pointer&lt; SAMRAI::hier::PatchHierarchy&lt; NDIM &gt; &gt; hierarchy, int level_number, double data_time, IBTK::LDataManager *l_data_manager)</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>IBAMR::IBStandardInitializer</name>
    <filename>class_i_b_a_m_r_1_1_i_b_standard_initializer.html</filename>
    <base>IBTK::LInitStrategy</base>
    <member kind="function">
      <type></type>
      <name>IBStandardInitializer</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_standard_initializer.html</anchorfile>
      <anchor>ad503c2a72ae0571f87a8e8d868f3bffe</anchor>
      <arglist>(const std::string &amp;object_name, SAMRAI::tbox::Pointer&lt; SAMRAI::tbox::Database &gt; input_db)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~IBStandardInitializer</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_standard_initializer.html</anchorfile>
      <anchor>a8a82ef71c4bb80fe788b8fa9a6a97367</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>registerLSiloDataWriter</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_standard_initializer.html</anchorfile>
      <anchor>ad7f149f596f5588f952391d9029566ae</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; IBTK::LSiloDataWriter &gt; silo_writer)</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>getLevelHasLagrangianData</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_standard_initializer.html</anchorfile>
      <anchor>add9b2b352ddfa78c0e7845e981d240ec</anchor>
      <arglist>(int level_number, bool can_be_refined) const </arglist>
    </member>
    <member kind="function">
      <type>unsigned int</type>
      <name>computeLocalNodeCountOnPatchLevel</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_standard_initializer.html</anchorfile>
      <anchor>a4c845b43b74cb36f9dc5a1227ce98821</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::hier::PatchHierarchy&lt; NDIM &gt; &gt; hierarchy, int level_number, double init_data_time, bool can_be_refined, bool initial_time)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>initializeStructureIndexingOnPatchLevel</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_standard_initializer.html</anchorfile>
      <anchor>a140796327416d7695d900ea4bfe709c0</anchor>
      <arglist>(std::map&lt; int, std::string &gt; &amp;strct_id_to_strct_name_map, std::map&lt; int, std::pair&lt; int, int &gt; &gt; &amp;strct_id_to_lag_idx_range_map, int level_number, double init_data_time, bool can_be_refined, bool initial_time, IBTK::LDataManager *l_data_manager)</arglist>
    </member>
    <member kind="function">
      <type>unsigned int</type>
      <name>initializeDataOnPatchLevel</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_standard_initializer.html</anchorfile>
      <anchor>a2623b3bd53c0a81d1fe77024d19bae4a</anchor>
      <arglist>(int lag_node_index_idx, unsigned int global_index_offset, unsigned int local_index_offset, SAMRAI::tbox::Pointer&lt; IBTK::LData &gt; X_data, SAMRAI::tbox::Pointer&lt; IBTK::LData &gt; U_data, SAMRAI::tbox::Pointer&lt; SAMRAI::hier::PatchHierarchy&lt; NDIM &gt; &gt; hierarchy, int level_number, double init_data_time, bool can_be_refined, bool initial_time, IBTK::LDataManager *l_data_manager)</arglist>
    </member>
    <member kind="function">
      <type>unsigned int</type>
      <name>initializeMassDataOnPatchLevel</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_standard_initializer.html</anchorfile>
      <anchor>afd73b60359d1ac7949fe9b088abe867a</anchor>
      <arglist>(unsigned int global_index_offset, unsigned int local_index_offset, SAMRAI::tbox::Pointer&lt; IBTK::LData &gt; M_data, SAMRAI::tbox::Pointer&lt; IBTK::LData &gt; K_data, SAMRAI::tbox::Pointer&lt; SAMRAI::hier::PatchHierarchy&lt; NDIM &gt; &gt; hierarchy, int level_number, double init_data_time, bool can_be_refined, bool initial_time, IBTK::LDataManager *l_data_manager)</arglist>
    </member>
    <member kind="function">
      <type>unsigned int</type>
      <name>initializeDirectorDataOnPatchLevel</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_standard_initializer.html</anchorfile>
      <anchor>a738dd9525290e6961874144b4dca3529</anchor>
      <arglist>(unsigned int global_index_offset, unsigned int local_index_offset, SAMRAI::tbox::Pointer&lt; IBTK::LData &gt; D_data, SAMRAI::tbox::Pointer&lt; SAMRAI::hier::PatchHierarchy&lt; NDIM &gt; &gt; hierarchy, int level_number, double init_data_time, bool can_be_refined, bool initial_time, IBTK::LDataManager *l_data_manager)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>tagCellsForInitialRefinement</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_standard_initializer.html</anchorfile>
      <anchor>a94188af9733365152857c2322fd0243a</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::hier::PatchHierarchy&lt; NDIM &gt; &gt; hierarchy, int level_number, double error_data_time, int tag_index)</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>IBAMR::IBStandardSourceGen</name>
    <filename>class_i_b_a_m_r_1_1_i_b_standard_source_gen.html</filename>
    <base>IBAMR::IBLagrangianSourceStrategy</base>
    <base>SAMRAI::tbox::Serializable</base>
    <member kind="function">
      <type></type>
      <name>IBStandardSourceGen</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_standard_source_gen.html</anchorfile>
      <anchor>ac93fa72f0504a82be7e5ea0beb86cb99</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~IBStandardSourceGen</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_standard_source_gen.html</anchorfile>
      <anchor>adafff45ab66fea1b422df722c7f70e22</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>std::vector&lt; double &gt; &amp;</type>
      <name>getSourceStrengths</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_standard_source_gen.html</anchorfile>
      <anchor>a6cd13c11f741b5b1883fd0e05b2e18fa</anchor>
      <arglist>(int ln)</arglist>
    </member>
    <member kind="function">
      <type>const std::vector&lt; double &gt; &amp;</type>
      <name>getSourceStrengths</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_standard_source_gen.html</anchorfile>
      <anchor>a0a1389f94759f7b9905f249dc55cab8a</anchor>
      <arglist>(int ln) const </arglist>
    </member>
    <member kind="function">
      <type>const std::vector&lt; double &gt; &amp;</type>
      <name>getSourcePressures</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_standard_source_gen.html</anchorfile>
      <anchor>a0ff27ec2753ec98efe27cfbc1f17a4b6</anchor>
      <arglist>(int ln) const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>initializeLevelData</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_standard_source_gen.html</anchorfile>
      <anchor>ad5314d332917ad4a0e112fef766c3eb1</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::hier::PatchHierarchy&lt; NDIM &gt; &gt; hierarchy, int level_number, double init_data_time, bool initial_time, IBTK::LDataManager *l_data_manager)</arglist>
    </member>
    <member kind="function">
      <type>unsigned int</type>
      <name>getNumSources</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_standard_source_gen.html</anchorfile>
      <anchor>a783644f05754051f18cbe0ee1a3f70a0</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::hier::PatchHierarchy&lt; NDIM &gt; &gt; hierarchy, int level_number, double data_time, IBTK::LDataManager *l_data_manager)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>getSourceLocations</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_standard_source_gen.html</anchorfile>
      <anchor>a045d09c30b1a4d284435aec18f651985</anchor>
      <arglist>(std::vector&lt; IBTK::Point &gt; &amp;X_src, std::vector&lt; double &gt; &amp;r_src, SAMRAI::tbox::Pointer&lt; IBTK::LData &gt; X_data, SAMRAI::tbox::Pointer&lt; SAMRAI::hier::PatchHierarchy&lt; NDIM &gt; &gt; hierarchy, int level_number, double data_time, IBTK::LDataManager *l_data_manager)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setSourcePressures</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_standard_source_gen.html</anchorfile>
      <anchor>ad28d18c8b15deb467f7b79e758b613ef</anchor>
      <arglist>(const std::vector&lt; double &gt; &amp;P_src, SAMRAI::tbox::Pointer&lt; SAMRAI::hier::PatchHierarchy&lt; NDIM &gt; &gt; hierarchy, int level_number, double data_time, IBTK::LDataManager *l_data_manager)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>computeSourceStrengths</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_standard_source_gen.html</anchorfile>
      <anchor>a204be665c4c899c194b9f5ea387f89af</anchor>
      <arglist>(std::vector&lt; double &gt; &amp;Q_src, SAMRAI::tbox::Pointer&lt; SAMRAI::hier::PatchHierarchy&lt; NDIM &gt; &gt; hierarchy, int level_number, double data_time, IBTK::LDataManager *l_data_manager)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>putToDatabase</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_standard_source_gen.html</anchorfile>
      <anchor>abbabe4b635f549320aa1cb5ce0b41cb3</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::tbox::Database &gt; db)</arglist>
    </member>
    <member kind="function" static="yes">
      <type>static bool</type>
      <name>getIsRegisteredWithStreamableManager</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_standard_source_gen.html</anchorfile>
      <anchor>a3f3a9be9af7e9abc4d09733eb7789f5b</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" static="yes">
      <type>static void</type>
      <name>setNumSources</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_standard_source_gen.html</anchorfile>
      <anchor>af800df86bfbe21049643a49249d4c069</anchor>
      <arglist>(int ln, unsigned int num_sources)</arglist>
    </member>
    <member kind="function" static="yes">
      <type>static unsigned int</type>
      <name>getNumSources</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_standard_source_gen.html</anchorfile>
      <anchor>a651797f69382145104c5e9b4696912f5</anchor>
      <arglist>(int ln)</arglist>
    </member>
    <member kind="function" static="yes">
      <type>static void</type>
      <name>setSourceNames</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_standard_source_gen.html</anchorfile>
      <anchor>aefdddd0ea848089a1c35d185e2fa5df5</anchor>
      <arglist>(int ln, const std::vector&lt; std::string &gt; &amp;names)</arglist>
    </member>
    <member kind="function" static="yes">
      <type>static const std::vector&lt; std::string &gt; &amp;</type>
      <name>getSourceNames</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_standard_source_gen.html</anchorfile>
      <anchor>a669d967907be9b651e6a6eed2b937922</anchor>
      <arglist>(int ln)</arglist>
    </member>
    <member kind="function" static="yes">
      <type>static void</type>
      <name>setSourceRadii</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_standard_source_gen.html</anchorfile>
      <anchor>a2c3675f4b95bf80af18bf2fdfa9885e0</anchor>
      <arglist>(int ln, const std::vector&lt; double &gt; &amp;radii)</arglist>
    </member>
    <member kind="function" static="yes">
      <type>static const std::vector&lt; double &gt; &amp;</type>
      <name>getSourceRadii</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_standard_source_gen.html</anchorfile>
      <anchor>acb2b48635d980a906c3d8705b8c912e3</anchor>
      <arglist>(int ln)</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>IBAMR::IBStrategy</name>
    <filename>class_i_b_a_m_r_1_1_i_b_strategy.html</filename>
    <base>StandardTagAndInitStrategy&lt; NDIM &gt;</base>
    <base>SAMRAI::tbox::Serializable</base>
  </compound>
  <compound kind="class">
    <name>IBAMR::IBStrategySet</name>
    <filename>class_i_b_a_m_r_1_1_i_b_strategy_set.html</filename>
    <base>IBAMR::IBStrategy</base>
    <member kind="function">
      <type></type>
      <name>IBStrategySet</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_strategy_set.html</anchorfile>
      <anchor>a09f980ec7dec1b7d2d348dc7dbf294a6</anchor>
      <arglist>(InputIterator first, InputIterator last)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~IBStrategySet</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_strategy_set.html</anchorfile>
      <anchor>afc7a9c1508b83e954da41f2afb322984</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>registerIBHierarchyIntegrator</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_strategy_set.html</anchorfile>
      <anchor>a6b7c526273c0d9e06e13afd4e7273c51</anchor>
      <arglist>(IBHierarchyIntegrator *ib_solver)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>registerEulerianVariables</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_strategy_set.html</anchorfile>
      <anchor>ace3a3d6f4f0919c0714bd011fc5868c6</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>registerEulerianCommunicationAlgorithms</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_strategy_set.html</anchorfile>
      <anchor>a5647a4c0188e420d1b1864f0d7a3b2ea</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>const SAMRAI::hier::IntVector&lt; NDIM &gt; &amp;</type>
      <name>getMinimumGhostCellWidth</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_strategy_set.html</anchorfile>
      <anchor>a1fbdd8e7209fa064d75c81eddd60e126</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setupTagBuffer</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_strategy_set.html</anchorfile>
      <anchor>a6c8af32b0328df65f4e5fe2b8e3109d9</anchor>
      <arglist>(SAMRAI::tbox::Array&lt; int &gt; &amp;tag_buffer, SAMRAI::tbox::Pointer&lt; SAMRAI::mesh::GriddingAlgorithm&lt; NDIM &gt; &gt; gridding_alg) const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>preprocessIntegrateData</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_strategy_set.html</anchorfile>
      <anchor>a9a02d8f849ce3f015631e58ab9a5eac3</anchor>
      <arglist>(double current_time, double new_time, int num_cycles)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>postprocessIntegrateData</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_strategy_set.html</anchorfile>
      <anchor>a52b7d605c2714a071ae686bcb4cce01d</anchor>
      <arglist>(double current_time, double new_time, int num_cycles)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>updateFixedLEOperators</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_strategy_set.html</anchorfile>
      <anchor>a2c547df850008917a0ff2a79d7d041fc</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>interpolateVelocity</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_strategy_set.html</anchorfile>
      <anchor>a24133f0568d718edb98d0d459c620408</anchor>
      <arglist>(int u_data_idx, const std::vector&lt; SAMRAI::tbox::Pointer&lt; SAMRAI::xfer::CoarsenSchedule&lt; NDIM &gt; &gt; &gt; &amp;u_synch_scheds, const std::vector&lt; SAMRAI::tbox::Pointer&lt; SAMRAI::xfer::RefineSchedule&lt; NDIM &gt; &gt; &gt; &amp;u_ghost_fill_scheds, double data_time)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>eulerStep</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_strategy_set.html</anchorfile>
      <anchor>ad9f0bbc75b2b16fcafcc96798e56cab9</anchor>
      <arglist>(double current_time, double new_time)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>midpointStep</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_strategy_set.html</anchorfile>
      <anchor>a83ba7e68cb6bf421c29be3c7ebff12a0</anchor>
      <arglist>(double current_time, double new_time)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>trapezoidalStep</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_strategy_set.html</anchorfile>
      <anchor>afe3fcf62c98a6ec7b50b71a05aafb849</anchor>
      <arglist>(double current_time, double new_time)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>computeLagrangianForce</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_strategy_set.html</anchorfile>
      <anchor>ae6ee96c60ab73ab5705b9fbeffc02cf0</anchor>
      <arglist>(double data_time)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>spreadForce</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_strategy_set.html</anchorfile>
      <anchor>aa3fd9fc6723be01193f3ee30298f9b95</anchor>
      <arglist>(int f_data_idx, IBTK::RobinPhysBdryPatchStrategy *f_phys_bdry_op, const std::vector&lt; SAMRAI::tbox::Pointer&lt; SAMRAI::xfer::RefineSchedule&lt; NDIM &gt; &gt; &gt; &amp;f_prolongation_scheds, double data_time)</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>hasFluidSources</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_strategy_set.html</anchorfile>
      <anchor>a02b4b3c12ffe3ae9c1d9057755b8ce9f</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>computeLagrangianFluidSource</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_strategy_set.html</anchorfile>
      <anchor>a01d9a7eea0908a0bd86470d5e4496529</anchor>
      <arglist>(double data_time)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>spreadFluidSource</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_strategy_set.html</anchorfile>
      <anchor>a6c79541f92d8152d6f65212e60f8446e</anchor>
      <arglist>(int q_data_idx, const std::vector&lt; SAMRAI::tbox::Pointer&lt; SAMRAI::xfer::RefineSchedule&lt; NDIM &gt; &gt; &gt; &amp;q_prolongation_scheds, double data_time)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>interpolatePressure</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_strategy_set.html</anchorfile>
      <anchor>a18e6828a7d1139e8598d60c833a38e8b</anchor>
      <arglist>(int p_data_idx, const std::vector&lt; SAMRAI::tbox::Pointer&lt; SAMRAI::xfer::CoarsenSchedule&lt; NDIM &gt; &gt; &gt; &amp;p_synch_scheds, const std::vector&lt; SAMRAI::tbox::Pointer&lt; SAMRAI::xfer::RefineSchedule&lt; NDIM &gt; &gt; &gt; &amp;p_ghost_fill_scheds, double data_time)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>preprocessSolveFluidEquations</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_strategy_set.html</anchorfile>
      <anchor>ad4381af0b88a83e5a6fae095f4e6b97e</anchor>
      <arglist>(double current_time, double new_time, int cycle_num)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>postprocessSolveFluidEquations</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_strategy_set.html</anchorfile>
      <anchor>a7077980c96448c6e7fc383d6cccb3432</anchor>
      <arglist>(double current_time, double new_time, int cycle_num)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>postprocessData</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_strategy_set.html</anchorfile>
      <anchor>a34b30690a23de603894c4aee33c2c653</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>initializePatchHierarchy</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_strategy_set.html</anchorfile>
      <anchor>aa06ffc9cfb46ed68561e88094ef75814</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::hier::PatchHierarchy&lt; NDIM &gt; &gt; hierarchy, SAMRAI::tbox::Pointer&lt; SAMRAI::mesh::GriddingAlgorithm&lt; NDIM &gt; &gt; gridding_alg, int u_data_idx, const std::vector&lt; SAMRAI::tbox::Pointer&lt; SAMRAI::xfer::CoarsenSchedule&lt; NDIM &gt; &gt; &gt; &amp;u_synch_scheds, const std::vector&lt; SAMRAI::tbox::Pointer&lt; SAMRAI::xfer::RefineSchedule&lt; NDIM &gt; &gt; &gt; &amp;u_ghost_fill_scheds, int integrator_step, double init_data_time, bool initial_time)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>registerLoadBalancer</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_strategy_set.html</anchorfile>
      <anchor>a60b8da00576780e1dce8c57020d827c1</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::mesh::LoadBalancer&lt; NDIM &gt; &gt; load_balancer, int workload_data_idx)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>updateWorkloadEstimates</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_strategy_set.html</anchorfile>
      <anchor>ad8812a205c848f64932cae44c75e48f2</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::hier::PatchHierarchy&lt; NDIM &gt; &gt; hierarchy, int workload_data_idx)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>beginDataRedistribution</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_strategy_set.html</anchorfile>
      <anchor>a57cc78c7019ceff02d79112b31420153</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::hier::PatchHierarchy&lt; NDIM &gt; &gt; hierarchy, SAMRAI::tbox::Pointer&lt; SAMRAI::mesh::GriddingAlgorithm&lt; NDIM &gt; &gt; gridding_alg)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>endDataRedistribution</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_strategy_set.html</anchorfile>
      <anchor>ac469e698a58768ccc95ac5c1619e1ae9</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::hier::PatchHierarchy&lt; NDIM &gt; &gt; hierarchy, SAMRAI::tbox::Pointer&lt; SAMRAI::mesh::GriddingAlgorithm&lt; NDIM &gt; &gt; gridding_alg)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>initializeLevelData</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_strategy_set.html</anchorfile>
      <anchor>a8a929a56ef86e980c2dd82aa5e0eacbb</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::hier::BasePatchHierarchy&lt; NDIM &gt; &gt; hierarchy, int level_number, double init_data_time, bool can_be_refined, bool initial_time, SAMRAI::tbox::Pointer&lt; SAMRAI::hier::BasePatchLevel&lt; NDIM &gt; &gt; old_level, bool allocate_data)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>resetHierarchyConfiguration</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_strategy_set.html</anchorfile>
      <anchor>a8039072023f06da8598a561e4d44196d</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::hier::BasePatchHierarchy&lt; NDIM &gt; &gt; hierarchy, int coarsest_level, int finest_level)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>applyGradientDetector</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_strategy_set.html</anchorfile>
      <anchor>ad53373eb6c5333b3b723b3d71471466a</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::hier::BasePatchHierarchy&lt; NDIM &gt; &gt; hierarchy, int level_number, double error_data_time, int tag_index, bool initial_time, bool uses_richardson_extrapolation_too)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>putToDatabase</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_strategy_set.html</anchorfile>
      <anchor>a09f1ad9d4f0b8ac3fa336f9b17c8b41e</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::tbox::Database &gt; db)</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>IBAMR::IBTargetPointForceSpec</name>
    <filename>class_i_b_a_m_r_1_1_i_b_target_point_force_spec.html</filename>
    <base>IBTK::Streamable</base>
    <member kind="function">
      <type></type>
      <name>IBTargetPointForceSpec</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_target_point_force_spec.html</anchorfile>
      <anchor>a90a2f597c658ebc7233b379114f506bb</anchor>
      <arglist>(int master_idx=-1, double kappa_target=0.0, double eta_target=0.0, const IBTK::Point &amp;X_target=IBTK::Point::Zero())</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~IBTargetPointForceSpec</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_target_point_force_spec.html</anchorfile>
      <anchor>a53e89240890ed9a373f018489d10e5d0</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>const int &amp;</type>
      <name>getMasterNodeIndex</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_target_point_force_spec.html</anchorfile>
      <anchor>af141105332ba1a91b5b85a94b0b0dad1</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>int &amp;</type>
      <name>getMasterNodeIndex</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_target_point_force_spec.html</anchorfile>
      <anchor>ac6a28ba73a8f0163e243adcb088ce171</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>const double &amp;</type>
      <name>getStiffness</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_target_point_force_spec.html</anchorfile>
      <anchor>a386e1b7edf36b5ad540bc35b389382d8</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>double &amp;</type>
      <name>getStiffness</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_target_point_force_spec.html</anchorfile>
      <anchor>a79fc1adcc560e9a586e2d04c74cf623e</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>const double &amp;</type>
      <name>getDamping</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_target_point_force_spec.html</anchorfile>
      <anchor>a7c3d73cfab8c4f53e1e7653606e93b36</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>double &amp;</type>
      <name>getDamping</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_target_point_force_spec.html</anchorfile>
      <anchor>a44546971a7f5479338d49ecc64863c93</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>const IBTK::Point &amp;</type>
      <name>getTargetPointPosition</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_target_point_force_spec.html</anchorfile>
      <anchor>a0b22b02a30165074f3047e3e94f4fe28</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>IBTK::Point &amp;</type>
      <name>getTargetPointPosition</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_target_point_force_spec.html</anchorfile>
      <anchor>a326c681b2904a7021605ea8abc9501eb</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>getStreamableClassID</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_target_point_force_spec.html</anchorfile>
      <anchor>a0b426b8e75042787b49ed8c4708e1811</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>size_t</type>
      <name>getDataStreamSize</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_target_point_force_spec.html</anchorfile>
      <anchor>aa226a4d1552eb6757396af968ee96628</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>packStream</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_target_point_force_spec.html</anchorfile>
      <anchor>a2d1ea5e6b0d492c24190f61d0d6520d9</anchor>
      <arglist>(SAMRAI::tbox::AbstractStream &amp;stream)</arglist>
    </member>
    <member kind="function" static="yes">
      <type>static void</type>
      <name>registerWithStreamableManager</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_target_point_force_spec.html</anchorfile>
      <anchor>a8773cbc3a1a5688b0e13305221818234</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" static="yes">
      <type>static bool</type>
      <name>getIsRegisteredWithStreamableManager</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_target_point_force_spec.html</anchorfile>
      <anchor>af97d5714893ad2d2be58071f83af9b09</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="variable" static="yes">
      <type>static int</type>
      <name>STREAMABLE_CLASS_ID</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_b_target_point_force_spec.html</anchorfile>
      <anchor>a5da993e930b54cee22a8dca91b701294</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>IBAMR::IMPInitializer</name>
    <filename>class_i_b_a_m_r_1_1_i_m_p_initializer.html</filename>
    <base>IBTK::LInitStrategy</base>
    <member kind="function">
      <type></type>
      <name>IMPInitializer</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_m_p_initializer.html</anchorfile>
      <anchor>a188ccf4112d6719af63406b78cdc10b7</anchor>
      <arglist>(const std::string &amp;object_name, SAMRAI::tbox::Pointer&lt; SAMRAI::tbox::Database &gt; input_db, SAMRAI::tbox::Pointer&lt; SAMRAI::hier::PatchHierarchy&lt; NDIM &gt; &gt; hierarchy, SAMRAI::tbox::Pointer&lt; SAMRAI::mesh::GriddingAlgorithm&lt; NDIM &gt; &gt; gridding_alg)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~IMPInitializer</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_m_p_initializer.html</anchorfile>
      <anchor>a50b8afb98a5a85154cbc07175d00ff8d</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>registerMesh</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_m_p_initializer.html</anchorfile>
      <anchor>af6a246df9418a7594fac3d1fbb052e76</anchor>
      <arglist>(libMesh::MeshBase *mesh, int level_number=-1)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>registerLSiloDataWriter</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_m_p_initializer.html</anchorfile>
      <anchor>a80fba00e10224f1b377e982072ea1329</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; IBTK::LSiloDataWriter &gt; silo_writer)</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>getLevelHasLagrangianData</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_m_p_initializer.html</anchorfile>
      <anchor>a2059b3422809bf653970f081866a0863</anchor>
      <arglist>(int level_number, bool can_be_refined) const </arglist>
    </member>
    <member kind="function">
      <type>unsigned int</type>
      <name>computeLocalNodeCountOnPatchLevel</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_m_p_initializer.html</anchorfile>
      <anchor>ac2723e7814fedc7a3a4ee4fc9ca5ed25</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::hier::PatchHierarchy&lt; NDIM &gt; &gt; hierarchy, int level_number, double init_data_time, bool can_be_refined, bool initial_time)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>initializeStructureIndexingOnPatchLevel</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_m_p_initializer.html</anchorfile>
      <anchor>ada0edcee02b1094f7d80d3b8ea2f8b36</anchor>
      <arglist>(std::map&lt; int, std::string &gt; &amp;strct_id_to_strct_name_map, std::map&lt; int, std::pair&lt; int, int &gt; &gt; &amp;strct_id_to_lag_idx_range_map, int level_number, double init_data_time, bool can_be_refined, bool initial_time, IBTK::LDataManager *l_data_manager)</arglist>
    </member>
    <member kind="function">
      <type>unsigned int</type>
      <name>initializeDataOnPatchLevel</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_m_p_initializer.html</anchorfile>
      <anchor>abe93e804bfcc182aa2b6862828e4c983</anchor>
      <arglist>(int lag_node_index_idx, unsigned int global_index_offset, unsigned int local_index_offset, SAMRAI::tbox::Pointer&lt; IBTK::LData &gt; X_data, SAMRAI::tbox::Pointer&lt; IBTK::LData &gt; U_data, SAMRAI::tbox::Pointer&lt; SAMRAI::hier::PatchHierarchy&lt; NDIM &gt; &gt; hierarchy, int level_number, double init_data_time, bool can_be_refined, bool initial_time, IBTK::LDataManager *l_data_manager)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>tagCellsForInitialRefinement</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_m_p_initializer.html</anchorfile>
      <anchor>af898f948a2217399d1b06188f92a78bc</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::hier::PatchHierarchy&lt; NDIM &gt; &gt; hierarchy, int level_number, double error_data_time, int tag_index)</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>IBAMR::IMPMethod</name>
    <filename>class_i_b_a_m_r_1_1_i_m_p_method.html</filename>
    <base>IBAMR::IBStrategy</base>
    <member kind="typedef">
      <type>void(*</type>
      <name>PK1StressFcnPtr</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_m_p_method.html</anchorfile>
      <anchor>ab7fa343275022fa9e0df5148664011ba</anchor>
      <arglist>)(libMesh::TensorValue&lt; double &gt; &amp;PP, const libMesh::TensorValue&lt; double &gt; &amp;FF, const libMesh::Point &amp;x, const libMesh::Point &amp;X, libMesh::subdomain_id_type subdomain_id, std::vector&lt; double &gt; &amp;internal_vars, double time, void *ctx)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>IMPMethod</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_m_p_method.html</anchorfile>
      <anchor>a4ef7584cfea6101ddda5715fc0b4243d</anchor>
      <arglist>(const std::string &amp;object_name, SAMRAI::tbox::Pointer&lt; SAMRAI::tbox::Database &gt; input_db, bool register_for_restart=true)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~IMPMethod</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_m_p_method.html</anchorfile>
      <anchor>abab8db5389bfbfb0478eb3bb4d4941cf</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>registerPK1StressTensorFunction</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_m_p_method.html</anchorfile>
      <anchor>a0c02ee3248159ab99c20accbe5c8ffdd</anchor>
      <arglist>(PK1StressFcnPtr PK1_stress_fcn, void *PK1_stress_fcn_ctx=NULL)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>registerLInitStrategy</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_m_p_method.html</anchorfile>
      <anchor>a6359b11f80c0ab7dc133ebd1c4d1d34e</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; IBTK::LInitStrategy &gt; l_initializer)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>freeLInitStrategy</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_m_p_method.html</anchorfile>
      <anchor>a400e75ae56d058020c533564ff9680b3</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>IBTK::LDataManager *</type>
      <name>getLDataManager</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_m_p_method.html</anchorfile>
      <anchor>a77aec1336ba857bbed2e4a4aa22dd2ef</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>registerLSiloDataWriter</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_m_p_method.html</anchorfile>
      <anchor>a6b9ed8ebc78a15a27e350bea63857fca</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; IBTK::LSiloDataWriter &gt; silo_writer)</arglist>
    </member>
    <member kind="function">
      <type>const SAMRAI::hier::IntVector&lt; NDIM &gt; &amp;</type>
      <name>getMinimumGhostCellWidth</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_m_p_method.html</anchorfile>
      <anchor>ae992d5e99ec74e47cc900bc272cc25f1</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setupTagBuffer</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_m_p_method.html</anchorfile>
      <anchor>a8fbbb762c4ddd19a9e22a85522e7cadf</anchor>
      <arglist>(SAMRAI::tbox::Array&lt; int &gt; &amp;tag_buffer, SAMRAI::tbox::Pointer&lt; SAMRAI::mesh::GriddingAlgorithm&lt; NDIM &gt; &gt; gridding_alg) const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>preprocessIntegrateData</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_m_p_method.html</anchorfile>
      <anchor>aebf29367ec644b211ec6eef08be4a039</anchor>
      <arglist>(double current_time, double new_time, int num_cycles)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>postprocessIntegrateData</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_m_p_method.html</anchorfile>
      <anchor>aea243d616fc1e3bc8a8cbba0bcf7dfe0</anchor>
      <arglist>(double current_time, double new_time, int num_cycles)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>interpolateVelocity</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_m_p_method.html</anchorfile>
      <anchor>aff1c67ab08c7fea9e92a6df312be05b2</anchor>
      <arglist>(int u_data_idx, const std::vector&lt; SAMRAI::tbox::Pointer&lt; SAMRAI::xfer::CoarsenSchedule&lt; NDIM &gt; &gt; &gt; &amp;u_synch_scheds, const std::vector&lt; SAMRAI::tbox::Pointer&lt; SAMRAI::xfer::RefineSchedule&lt; NDIM &gt; &gt; &gt; &amp;u_ghost_fill_scheds, double data_time)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>eulerStep</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_m_p_method.html</anchorfile>
      <anchor>a2bcd589bf28943cd8f0f2c84f38ce95a</anchor>
      <arglist>(double current_time, double new_time)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>midpointStep</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_m_p_method.html</anchorfile>
      <anchor>a53ecccd06db1780401a5b0fbff3a58b1</anchor>
      <arglist>(double current_time, double new_time)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>trapezoidalStep</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_m_p_method.html</anchorfile>
      <anchor>a622546b97262b54a74c7593e067f0638</anchor>
      <arglist>(double current_time, double new_time)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>computeLagrangianForce</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_m_p_method.html</anchorfile>
      <anchor>a9a97e28cb0d853d52a6903eed4f10732</anchor>
      <arglist>(double data_time)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>spreadForce</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_m_p_method.html</anchorfile>
      <anchor>ab6db2b78d7d8e85632f12f006cc271fe</anchor>
      <arglist>(int f_data_idx, IBTK::RobinPhysBdryPatchStrategy *f_phys_bdry_op, const std::vector&lt; SAMRAI::tbox::Pointer&lt; SAMRAI::xfer::RefineSchedule&lt; NDIM &gt; &gt; &gt; &amp;f_prolongation_scheds, double data_time)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>initializePatchHierarchy</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_m_p_method.html</anchorfile>
      <anchor>a52de0f65c5f9af9ec17da557d4409c91</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::hier::PatchHierarchy&lt; NDIM &gt; &gt; hierarchy, SAMRAI::tbox::Pointer&lt; SAMRAI::mesh::GriddingAlgorithm&lt; NDIM &gt; &gt; gridding_alg, int u_data_idx, const std::vector&lt; SAMRAI::tbox::Pointer&lt; SAMRAI::xfer::CoarsenSchedule&lt; NDIM &gt; &gt; &gt; &amp;u_synch_scheds, const std::vector&lt; SAMRAI::tbox::Pointer&lt; SAMRAI::xfer::RefineSchedule&lt; NDIM &gt; &gt; &gt; &amp;u_ghost_fill_scheds, int integrator_step, double init_data_time, bool initial_time)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>registerLoadBalancer</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_m_p_method.html</anchorfile>
      <anchor>a7db6ca79f3461a246732c01b05bd1288</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::mesh::LoadBalancer&lt; NDIM &gt; &gt; load_balancer, int workload_data_idx)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>updateWorkloadEstimates</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_m_p_method.html</anchorfile>
      <anchor>a5e3463db0456d9bbebbd655973f44c76</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::hier::PatchHierarchy&lt; NDIM &gt; &gt; hierarchy, int workload_data_idx)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>beginDataRedistribution</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_m_p_method.html</anchorfile>
      <anchor>a9f0702b03fddb9cdc7afc407cd2a034a</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::hier::PatchHierarchy&lt; NDIM &gt; &gt; hierarchy, SAMRAI::tbox::Pointer&lt; SAMRAI::mesh::GriddingAlgorithm&lt; NDIM &gt; &gt; gridding_alg)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>endDataRedistribution</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_m_p_method.html</anchorfile>
      <anchor>a8802a721daefc467a83abf900a6f4079</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::hier::PatchHierarchy&lt; NDIM &gt; &gt; hierarchy, SAMRAI::tbox::Pointer&lt; SAMRAI::mesh::GriddingAlgorithm&lt; NDIM &gt; &gt; gridding_alg)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>initializeLevelData</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_m_p_method.html</anchorfile>
      <anchor>a318267d6e2411b8235417890e4a5a14f</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::hier::BasePatchHierarchy&lt; NDIM &gt; &gt; hierarchy, int level_number, double init_data_time, bool can_be_refined, bool initial_time, SAMRAI::tbox::Pointer&lt; SAMRAI::hier::BasePatchLevel&lt; NDIM &gt; &gt; old_level, bool allocate_data)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>resetHierarchyConfiguration</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_m_p_method.html</anchorfile>
      <anchor>af9927e64da21a1cfa652b9ab21b742ae</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::hier::BasePatchHierarchy&lt; NDIM &gt; &gt; hierarchy, int coarsest_level, int finest_level)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>applyGradientDetector</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_m_p_method.html</anchorfile>
      <anchor>a9a1702da549874b2b62b2dc6887ec91d</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::hier::BasePatchHierarchy&lt; NDIM &gt; &gt; hierarchy, int level_number, double error_data_time, int tag_index, bool initial_time, bool uses_richardson_extrapolation_too)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>putToDatabase</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_m_p_method.html</anchorfile>
      <anchor>aea8afe044bcdc631798942295bb54650</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::tbox::Database &gt; db)</arglist>
    </member>
    <member kind="function" protection="protected">
      <type>void</type>
      <name>getPositionData</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_m_p_method.html</anchorfile>
      <anchor>a6c49f42969d6442e3928234021974d79</anchor>
      <arglist>(std::vector&lt; SAMRAI::tbox::Pointer&lt; IBTK::LData &gt; &gt; **X_data, bool **X_needs_ghost_fill, double data_time)</arglist>
    </member>
    <member kind="function" protection="protected">
      <type>void</type>
      <name>getVelocityData</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_m_p_method.html</anchorfile>
      <anchor>aacaa8fa99c89ad677dfbd6ce11d63ade</anchor>
      <arglist>(std::vector&lt; SAMRAI::tbox::Pointer&lt; IBTK::LData &gt; &gt; **U_data, std::vector&lt; SAMRAI::tbox::Pointer&lt; IBTK::LData &gt; &gt; **Grad_U_data, double data_time)</arglist>
    </member>
    <member kind="function" protection="protected">
      <type>void</type>
      <name>getDeformationGradientData</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_m_p_method.html</anchorfile>
      <anchor>a0a4a67cba3ba5c26590c2fd4ca59081b</anchor>
      <arglist>(std::vector&lt; SAMRAI::tbox::Pointer&lt; IBTK::LData &gt; &gt; **F_data, double data_time)</arglist>
    </member>
    <member kind="function" protection="protected">
      <type>void</type>
      <name>reinitMidpointData</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_m_p_method.html</anchorfile>
      <anchor>a593fabb6a03ec40c4f3559c8b8c54234</anchor>
      <arglist>(const std::vector&lt; SAMRAI::tbox::Pointer&lt; IBTK::LData &gt; &gt; &amp;current_data, const std::vector&lt; SAMRAI::tbox::Pointer&lt; IBTK::LData &gt; &gt; &amp;new_data, const std::vector&lt; SAMRAI::tbox::Pointer&lt; IBTK::LData &gt; &gt; &amp;half_data)</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>IBAMR::INSCollocatedCenteredConvectiveOperator</name>
    <filename>class_i_b_a_m_r_1_1_i_n_s_collocated_centered_convective_operator.html</filename>
    <base>IBAMR::ConvectiveOperator</base>
    <member kind="function">
      <type></type>
      <name>INSCollocatedCenteredConvectiveOperator</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_n_s_collocated_centered_convective_operator.html</anchorfile>
      <anchor>a43027648c10d0929b6e2064c3759809c</anchor>
      <arglist>(const std::string &amp;object_name, SAMRAI::tbox::Pointer&lt; SAMRAI::tbox::Database &gt; input_db, ConvectiveDifferencingType difference_form, const std::vector&lt; SAMRAI::solv::RobinBcCoefStrategy&lt; NDIM &gt; * &gt; &amp;bc_coefs)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~INSCollocatedCenteredConvectiveOperator</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_n_s_collocated_centered_convective_operator.html</anchorfile>
      <anchor>a6823f102e9fb860f6cd6205bcfeaaf37</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>applyConvectiveOperator</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_n_s_collocated_centered_convective_operator.html</anchorfile>
      <anchor>a6401b22116918c32a17d12fd0463d792</anchor>
      <arglist>(int U_idx, int N_idx)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>initializeOperatorState</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_n_s_collocated_centered_convective_operator.html</anchorfile>
      <anchor>ab9a983921ec1809540aab5a58585e4c5</anchor>
      <arglist>(const SAMRAI::solv::SAMRAIVectorReal&lt; NDIM, double &gt; &amp;in, const SAMRAI::solv::SAMRAIVectorReal&lt; NDIM, double &gt; &amp;out)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>deallocateOperatorState</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_n_s_collocated_centered_convective_operator.html</anchorfile>
      <anchor>a41b8b8e13d3eb9e65eca0e697860fd19</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" static="yes">
      <type>static SAMRAI::tbox::Pointer&lt; ConvectiveOperator &gt;</type>
      <name>allocate_operator</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_n_s_collocated_centered_convective_operator.html</anchorfile>
      <anchor>a7ce20912d57f623b4f38d546a9a13b7b</anchor>
      <arglist>(const std::string &amp;object_name, SAMRAI::tbox::Pointer&lt; SAMRAI::tbox::Database &gt; input_db, ConvectiveDifferencingType difference_form, const std::vector&lt; SAMRAI::solv::RobinBcCoefStrategy&lt; NDIM &gt; * &gt; &amp;bc_coefs)</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>IBAMR::INSCollocatedConvectiveOperatorManager</name>
    <filename>class_i_b_a_m_r_1_1_i_n_s_collocated_convective_operator_manager.html</filename>
    <member kind="typedef">
      <type>SAMRAI::tbox::Pointer&lt; ConvectiveOperator &gt;(*</type>
      <name>OperatorMaker</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_n_s_collocated_convective_operator_manager.html</anchorfile>
      <anchor>a3fc3ded77b8fd7de979eb7b0cc030161</anchor>
      <arglist>)(const std::string &amp;operator_object_name, SAMRAI::tbox::Pointer&lt; SAMRAI::tbox::Database &gt; input_db, ConvectiveDifferencingType difference_form, const std::vector&lt; SAMRAI::solv::RobinBcCoefStrategy&lt; NDIM &gt; * &gt; &amp;bc_coefs)</arglist>
    </member>
    <member kind="function">
      <type>SAMRAI::tbox::Pointer&lt; ConvectiveOperator &gt;</type>
      <name>allocateOperator</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_n_s_collocated_convective_operator_manager.html</anchorfile>
      <anchor>aeb013dadb4561cc9c9acde052dec9ad7</anchor>
      <arglist>(const std::string &amp;operator_type, const std::string &amp;operator_object_name, SAMRAI::tbox::Pointer&lt; SAMRAI::tbox::Database &gt; input_db, ConvectiveDifferencingType difference_form, const std::vector&lt; SAMRAI::solv::RobinBcCoefStrategy&lt; NDIM &gt; * &gt; &amp;bc_coefs) const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>registerOperatorFactoryFunction</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_n_s_collocated_convective_operator_manager.html</anchorfile>
      <anchor>a9278c4f61e0e201158d4d16ee51b7ad1</anchor>
      <arglist>(const std::string &amp;operator_type, OperatorMaker operator_maker)</arglist>
    </member>
    <member kind="function" static="yes">
      <type>static INSCollocatedConvectiveOperatorManager *</type>
      <name>getManager</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_n_s_collocated_convective_operator_manager.html</anchorfile>
      <anchor>a57c391f099e3cdb6ec24727430e410d4</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" static="yes">
      <type>static void</type>
      <name>freeManager</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_n_s_collocated_convective_operator_manager.html</anchorfile>
      <anchor>a9ec65aca7ffedad9a1232d3474b32d19</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="variable" static="yes">
      <type>static const std::string</type>
      <name>DEFAULT</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_n_s_collocated_convective_operator_manager.html</anchorfile>
      <anchor>a325c6ffb7856aeee28eb6a9a70cc6f73</anchor>
      <arglist></arglist>
    </member>
    <member kind="function" protection="protected">
      <type></type>
      <name>INSCollocatedConvectiveOperatorManager</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_n_s_collocated_convective_operator_manager.html</anchorfile>
      <anchor>a04748098eb0ea802e8e7f34687df3f79</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" protection="protected">
      <type></type>
      <name>~INSCollocatedConvectiveOperatorManager</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_n_s_collocated_convective_operator_manager.html</anchorfile>
      <anchor>a1b9e85774cd9ba7d40649b97ce3befe2</anchor>
      <arglist>()</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>IBAMR::INSCollocatedHierarchyIntegrator</name>
    <filename>class_i_b_a_m_r_1_1_i_n_s_collocated_hierarchy_integrator.html</filename>
    <base>IBAMR::INSHierarchyIntegrator</base>
    <member kind="function">
      <type></type>
      <name>INSCollocatedHierarchyIntegrator</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_n_s_collocated_hierarchy_integrator.html</anchorfile>
      <anchor>a0e3a40200e12fff0c6ffa18e3690baaa</anchor>
      <arglist>(const std::string &amp;object_name, SAMRAI::tbox::Pointer&lt; SAMRAI::tbox::Database &gt; input_db, bool register_for_restart=true)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~INSCollocatedHierarchyIntegrator</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_n_s_collocated_hierarchy_integrator.html</anchorfile>
      <anchor>a684a003d81335996e9d2d4c8d411d8eb</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>SAMRAI::tbox::Pointer&lt; ConvectiveOperator &gt;</type>
      <name>getConvectiveOperator</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_n_s_collocated_hierarchy_integrator.html</anchorfile>
      <anchor>afbd3103cb5aac8f9e3344caea8f43c12</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>SAMRAI::tbox::Pointer&lt; IBTK::PoissonSolver &gt;</type>
      <name>getVelocitySubdomainSolver</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_n_s_collocated_hierarchy_integrator.html</anchorfile>
      <anchor>aad43a1d1731f81b2c29178fb59a85342</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>SAMRAI::tbox::Pointer&lt; IBTK::PoissonSolver &gt;</type>
      <name>getPressureSubdomainSolver</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_n_s_collocated_hierarchy_integrator.html</anchorfile>
      <anchor>a7ccc4882966c4615eb369d2a6638f26e</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>initializeHierarchyIntegrator</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_n_s_collocated_hierarchy_integrator.html</anchorfile>
      <anchor>a955179ef975d1fc9526709348f8b2ee3</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::hier::PatchHierarchy&lt; NDIM &gt; &gt; hierarchy, SAMRAI::tbox::Pointer&lt; SAMRAI::mesh::GriddingAlgorithm&lt; NDIM &gt; &gt; gridding_alg)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>initializePatchHierarchy</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_n_s_collocated_hierarchy_integrator.html</anchorfile>
      <anchor>a809484577f7843097f10a981f720eb23</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::hier::PatchHierarchy&lt; NDIM &gt; &gt; hierarchy, SAMRAI::tbox::Pointer&lt; SAMRAI::mesh::GriddingAlgorithm&lt; NDIM &gt; &gt; gridding_alg)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>preprocessIntegrateHierarchy</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_n_s_collocated_hierarchy_integrator.html</anchorfile>
      <anchor>a739a69560e9e194610ba47f6c098338e</anchor>
      <arglist>(double current_time, double new_time, int num_cycles=1)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>integrateHierarchy</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_n_s_collocated_hierarchy_integrator.html</anchorfile>
      <anchor>a9995f3bf771e37d62e35238517d8489c</anchor>
      <arglist>(double current_time, double new_time, int cycle_num=0)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>postprocessIntegrateHierarchy</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_n_s_collocated_hierarchy_integrator.html</anchorfile>
      <anchor>ab7c1bb352b13da30765147580fd8a5fa</anchor>
      <arglist>(double current_time, double new_time, bool skip_synchronize_new_state_data, int num_cycles=1)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>regridHierarchy</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_n_s_collocated_hierarchy_integrator.html</anchorfile>
      <anchor>abe7cf638d770b069122d416ccfa232fe</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~INSHierarchyIntegrator</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_n_s_hierarchy_integrator.html</anchorfile>
      <anchor>aee3b8c472de36a645b38742b2dfaa34c</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setViscousTimeSteppingType</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_n_s_hierarchy_integrator.html</anchorfile>
      <anchor>a0117fab736464729e16f9e674191de7f</anchor>
      <arglist>(TimeSteppingType viscous_time_stepping_type)</arglist>
    </member>
    <member kind="function">
      <type>TimeSteppingType</type>
      <name>getViscousTimeSteppingType</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_n_s_hierarchy_integrator.html</anchorfile>
      <anchor>ab6964f88a05d76f92a501f0cc3bb0bb8</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setConvectiveTimeSteppingType</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_n_s_hierarchy_integrator.html</anchorfile>
      <anchor>a87fe194e158ba1111bd8dda21d1f5a9c</anchor>
      <arglist>(TimeSteppingType convective_time_stepping_type)</arglist>
    </member>
    <member kind="function">
      <type>TimeSteppingType</type>
      <name>getConvectiveTimeSteppingType</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_n_s_hierarchy_integrator.html</anchorfile>
      <anchor>a4959580e78a73ce9ee1e822444ce4449</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setInitialConvectiveTimeSteppingType</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_n_s_hierarchy_integrator.html</anchorfile>
      <anchor>abb0a8fce484cf90109633e6d61140374</anchor>
      <arglist>(TimeSteppingType init_convective_time_stepping_type) const </arglist>
    </member>
    <member kind="function">
      <type>TimeSteppingType</type>
      <name>getInitialConvectiveTimeSteppingType</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_n_s_hierarchy_integrator.html</anchorfile>
      <anchor>a54b865f40e7790780f473fa2b4118bc9</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>registerAdvDiffHierarchyIntegrator</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_n_s_hierarchy_integrator.html</anchorfile>
      <anchor>a1f029dcf5956a842383a579ecc2decd0</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; AdvDiffHierarchyIntegrator &gt; adv_diff_hier_integrator)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setStokesSpecifications</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_n_s_hierarchy_integrator.html</anchorfile>
      <anchor>a4575ad2ae4bf529b04d9b8d01e3427ba</anchor>
      <arglist>(StokesSpecifications problem_coefs)</arglist>
    </member>
    <member kind="function">
      <type>const StokesSpecifications *</type>
      <name>getStokesSpecifications</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_n_s_hierarchy_integrator.html</anchorfile>
      <anchor>a14fe01914085e26d15e8bc08cd7c16d3</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>registerPhysicalBoundaryConditions</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_n_s_hierarchy_integrator.html</anchorfile>
      <anchor>a3751b764c51e96a286d234b3d0e9023e</anchor>
      <arglist>(const std::vector&lt; SAMRAI::solv::RobinBcCoefStrategy&lt; NDIM &gt; * &gt; &amp;bc_coefs)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual const std::vector&lt; SAMRAI::solv::RobinBcCoefStrategy&lt; NDIM &gt; * &gt; &amp;</type>
      <name>getVelocityBoundaryConditions</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_n_s_hierarchy_integrator.html</anchorfile>
      <anchor>ab3a58da2c5b934a54f69ab220988311c</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual SAMRAI::solv::RobinBcCoefStrategy&lt; NDIM &gt; *</type>
      <name>getPressureBoundaryConditions</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_n_s_hierarchy_integrator.html</anchorfile>
      <anchor>a04848045a08cb1fcb8e3795da7002fad</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>registerVelocityInitialConditions</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_n_s_hierarchy_integrator.html</anchorfile>
      <anchor>af74c25be347bbb133a9f37bbfe44de1b</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; IBTK::CartGridFunction &gt; U_init)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>registerPressureInitialConditions</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_n_s_hierarchy_integrator.html</anchorfile>
      <anchor>a08c932b300ca7946dcd8bca37f4d18c9</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; IBTK::CartGridFunction &gt; P_init)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>registerBodyForceFunction</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_n_s_hierarchy_integrator.html</anchorfile>
      <anchor>af7d95d823afa9aecca1662c793ce54aa</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; IBTK::CartGridFunction &gt; F_fcn)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>registerFluidSourceFunction</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_n_s_hierarchy_integrator.html</anchorfile>
      <anchor>a3519c1d326386ea21466904aeaf18925</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; IBTK::CartGridFunction &gt; Q_fcn)</arglist>
    </member>
    <member kind="function">
      <type>SAMRAI::tbox::Pointer&lt; SAMRAI::hier::Variable&lt; NDIM &gt; &gt;</type>
      <name>getVelocityVariable</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_n_s_hierarchy_integrator.html</anchorfile>
      <anchor>acc99d272b262c2d526dbaaa0fd0caf31</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>SAMRAI::tbox::Pointer&lt; SAMRAI::hier::Variable&lt; NDIM &gt; &gt;</type>
      <name>getPressureVariable</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_n_s_hierarchy_integrator.html</anchorfile>
      <anchor>a56988073b107f99fd0021d48f2ffd617</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>SAMRAI::tbox::Pointer&lt; SAMRAI::hier::Variable&lt; NDIM &gt; &gt;</type>
      <name>getBodyForceVariable</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_n_s_hierarchy_integrator.html</anchorfile>
      <anchor>a7a3884f26f7b95c029bafde430f4af3c</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>SAMRAI::tbox::Pointer&lt; SAMRAI::hier::Variable&lt; NDIM &gt; &gt;</type>
      <name>getFluidSourceVariable</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_n_s_hierarchy_integrator.html</anchorfile>
      <anchor>a0aa6680357f6951e8cde8d12fd00dd21</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::FaceVariable&lt; NDIM, double &gt; &gt;</type>
      <name>getAdvectionVelocityVariable</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_n_s_hierarchy_integrator.html</anchorfile>
      <anchor>ac27dd505708d24b0ed678b8b5fa5081b</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>std::vector&lt; SAMRAI::solv::RobinBcCoefStrategy&lt; NDIM &gt; * &gt;</type>
      <name>getIntermediateVelocityBoundaryConditions</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_n_s_hierarchy_integrator.html</anchorfile>
      <anchor>a210293125183f4d79735ca7a74ebd786</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>SAMRAI::solv::RobinBcCoefStrategy&lt; NDIM &gt; *</type>
      <name>getProjectionBoundaryConditions</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_n_s_hierarchy_integrator.html</anchorfile>
      <anchor>a5d5d3ca230e17ccc4c7fcb5342f808cd</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>registerMassDensityVariable</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_n_s_hierarchy_integrator.html</anchorfile>
      <anchor>a9f3b397d1b06963287c97d9216c7589e</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::hier::Variable&lt; NDIM &gt; &gt; rho_var)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setMassDensityFunction</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_n_s_hierarchy_integrator.html</anchorfile>
      <anchor>a6ce05ca032447936b33cb219b85220cf</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; IBTK::CartGridFunction &gt; rho_fcn)</arglist>
    </member>
    <member kind="function">
      <type>SAMRAI::tbox::Pointer&lt; IBTK::CartGridFunction &gt;</type>
      <name>getMassDensityFunction</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_n_s_hierarchy_integrator.html</anchorfile>
      <anchor>a3734e7fdc14cb92ba5464b0f85cd0b13</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setConvectiveOperatorType</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_n_s_hierarchy_integrator.html</anchorfile>
      <anchor>aa22f4b5dbe32a30dd8cf191ff5a7399f</anchor>
      <arglist>(const std::string &amp;op_type)</arglist>
    </member>
    <member kind="function">
      <type>const std::string &amp;</type>
      <name>getConvectiveOperatorType</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_n_s_hierarchy_integrator.html</anchorfile>
      <anchor>af43adc13e5b65b93437be5e5b48372f4</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setConvectiveDifferencingType</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_n_s_hierarchy_integrator.html</anchorfile>
      <anchor>a243c1af224d4f2934f971267407fb209</anchor>
      <arglist>(ConvectiveDifferencingType difference_form)</arglist>
    </member>
    <member kind="function">
      <type>ConvectiveDifferencingType</type>
      <name>getConvectiveDifferencingType</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_n_s_hierarchy_integrator.html</anchorfile>
      <anchor>a8db46f081ea59988192f3435e32641ae</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setCreepingFlow</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_n_s_hierarchy_integrator.html</anchorfile>
      <anchor>ade3e2dffa2d179e033fa52242bdf7e2e</anchor>
      <arglist>(bool creeping_flow)</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>getCreepingFlow</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_n_s_hierarchy_integrator.html</anchorfile>
      <anchor>a5e31031d4d078427e66d532114451a28</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setConvectiveOperator</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_n_s_hierarchy_integrator.html</anchorfile>
      <anchor>ac133a7e601e622b15d15160d7b27a569</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; ConvectiveOperator &gt; convective_op)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setVelocitySubdomainSolver</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_n_s_hierarchy_integrator.html</anchorfile>
      <anchor>aa9864fa11242c9a043d9e3f4383bdcbb</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; IBTK::PoissonSolver &gt; velocity_solver)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setPressureSubdomainSolver</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_n_s_hierarchy_integrator.html</anchorfile>
      <anchor>aa2c88006d9237e43666153b3c6db408f</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; IBTK::PoissonSolver &gt; pressure_solver)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>getNumberOfCycles</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_n_s_hierarchy_integrator.html</anchorfile>
      <anchor>a4cf6c5ec22c689026de6f0c411befa0d</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function" protection="protected">
      <type>double</type>
      <name>getStableTimestep</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_n_s_collocated_hierarchy_integrator.html</anchorfile>
      <anchor>aa6023ffa202f63bc8328de6705d9ecfe</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::hier::Patch&lt; NDIM &gt; &gt; patch) const </arglist>
    </member>
    <member kind="function" protection="protected">
      <type>void</type>
      <name>initializeLevelDataSpecialized</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_n_s_collocated_hierarchy_integrator.html</anchorfile>
      <anchor>a85982ad3e5e6eeae528a3d699d77f213</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::hier::BasePatchHierarchy&lt; NDIM &gt; &gt; hierarchy, int level_number, double init_data_time, bool can_be_refined, bool initial_time, SAMRAI::tbox::Pointer&lt; SAMRAI::hier::BasePatchLevel&lt; NDIM &gt; &gt; old_level, bool allocate_data)</arglist>
    </member>
    <member kind="function" protection="protected">
      <type>void</type>
      <name>resetHierarchyConfigurationSpecialized</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_n_s_collocated_hierarchy_integrator.html</anchorfile>
      <anchor>af57f9e7e45e79df3bd8565302561df73</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::hier::BasePatchHierarchy&lt; NDIM &gt; &gt; hierarchy, int coarsest_level, int finest_level)</arglist>
    </member>
    <member kind="function" protection="protected">
      <type>void</type>
      <name>applyGradientDetectorSpecialized</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_n_s_collocated_hierarchy_integrator.html</anchorfile>
      <anchor>ad8fac75131cdae16ae4f3e0209cbc4b6</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::hier::BasePatchHierarchy&lt; NDIM &gt; &gt; hierarchy, int level_number, double error_data_time, int tag_index, bool initial_time, bool uses_richardson_extrapolation_too)</arglist>
    </member>
    <member kind="function" protection="protected">
      <type>void</type>
      <name>setupPlotDataSpecialized</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_n_s_collocated_hierarchy_integrator.html</anchorfile>
      <anchor>a9d5703a1be086e65e918a3d277b0804d</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" protection="protected">
      <type></type>
      <name>INSHierarchyIntegrator</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_n_s_hierarchy_integrator.html</anchorfile>
      <anchor>aeae6eb9f7af3ccf03c3d8e9e788f1b12</anchor>
      <arglist>(const std::string &amp;object_name, SAMRAI::tbox::Pointer&lt; SAMRAI::tbox::Database &gt; input_db, SAMRAI::tbox::Pointer&lt; SAMRAI::hier::Variable&lt; NDIM &gt; &gt; U_var, SAMRAI::tbox::Pointer&lt; SAMRAI::hier::Variable&lt; NDIM &gt; &gt; P_var, SAMRAI::tbox::Pointer&lt; SAMRAI::hier::Variable&lt; NDIM &gt; &gt; F_var, SAMRAI::tbox::Pointer&lt; SAMRAI::hier::Variable&lt; NDIM &gt; &gt; Q_var, bool register_for_restart)</arglist>
    </member>
    <member kind="function" protection="protected">
      <type>double</type>
      <name>getMaximumTimeStepSizeSpecialized</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_n_s_hierarchy_integrator.html</anchorfile>
      <anchor>a78b448dc32ee20dd85dfe82e26d7ebf3</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" protection="protected">
      <type>double</type>
      <name>getStableTimestep</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_n_s_hierarchy_integrator.html</anchorfile>
      <anchor>aadab46eac1f9da13b441dd31277953e6</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::hier::PatchLevel&lt; NDIM &gt; &gt; level) const </arglist>
    </member>
    <member kind="function" protection="protected">
      <type>void</type>
      <name>putToDatabaseSpecialized</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_n_s_hierarchy_integrator.html</anchorfile>
      <anchor>a89829f9e438e9f98f2bda75c5d03cbc6</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::tbox::Database &gt; db)</arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>TimeSteppingType</type>
      <name>d_viscous_time_stepping_type</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_n_s_hierarchy_integrator.html</anchorfile>
      <anchor>aa77e68062fcd5ec592a2304ebaf0376f</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>TimeSteppingType</type>
      <name>d_convective_time_stepping_type</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_n_s_hierarchy_integrator.html</anchorfile>
      <anchor>ab87eaa443308b1829c48546b22d56a63</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>TimeSteppingType</type>
      <name>d_init_convective_time_stepping_type</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_n_s_hierarchy_integrator.html</anchorfile>
      <anchor>a7720811347ef2fbfb15dcee4a50f4b07</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>StokesSpecifications</type>
      <name>d_problem_coefs</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_n_s_hierarchy_integrator.html</anchorfile>
      <anchor>aa428d2ebcf601eccc0068aefb9a7c175</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>double</type>
      <name>d_cfl_max</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_n_s_hierarchy_integrator.html</anchorfile>
      <anchor>a5210a458e89f6d386fffa7f99f9432fe</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>bool</type>
      <name>d_using_vorticity_tagging</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_n_s_hierarchy_integrator.html</anchorfile>
      <anchor>a5b1c616621428979719b88c88ba6446c</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>bool</type>
      <name>d_normalize_pressure</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_n_s_hierarchy_integrator.html</anchorfile>
      <anchor>a4e70b8ced04e27b33e40a4d489767a1a</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>bool</type>
      <name>d_normalize_velocity</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_n_s_hierarchy_integrator.html</anchorfile>
      <anchor>af8654e0d7128c474d6adc04597a9f274</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>bool</type>
      <name>d_creeping_flow</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_n_s_hierarchy_integrator.html</anchorfile>
      <anchor>a8deab32b94a767664f670479ec7dfd66</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>double</type>
      <name>d_regrid_max_div_growth_factor</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_n_s_hierarchy_integrator.html</anchorfile>
      <anchor>ac708c7655be7ebd021842b38406495a7</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>double</type>
      <name>d_U_scale</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_n_s_hierarchy_integrator.html</anchorfile>
      <anchor>ad19b2e748227ca58c5e024643b163484</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>SAMRAI::tbox::Pointer&lt; SAMRAI::hier::Variable&lt; NDIM &gt; &gt;</type>
      <name>d_U_var</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_n_s_hierarchy_integrator.html</anchorfile>
      <anchor>ac82eea36d04c6b16da6a95be71126cd8</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>SAMRAI::tbox::Pointer&lt; IBTK::CartGridFunction &gt;</type>
      <name>d_U_init</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_n_s_hierarchy_integrator.html</anchorfile>
      <anchor>acd3a4571780df4a4cd2964d60a87e417</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>int</type>
      <name>d_coarsest_reset_ln</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_n_s_hierarchy_integrator.html</anchorfile>
      <anchor>adb4b0289e8699f18e581be24c096f3dc</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>IBAMR::INSCollocatedPPMConvectiveOperator</name>
    <filename>class_i_b_a_m_r_1_1_i_n_s_collocated_p_p_m_convective_operator.html</filename>
    <base>IBAMR::ConvectiveOperator</base>
    <member kind="function">
      <type></type>
      <name>INSCollocatedPPMConvectiveOperator</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_n_s_collocated_p_p_m_convective_operator.html</anchorfile>
      <anchor>ab58b66ca1e6823671adf8e954fd6283d</anchor>
      <arglist>(const std::string &amp;object_name, SAMRAI::tbox::Pointer&lt; SAMRAI::tbox::Database &gt; input_db, ConvectiveDifferencingType difference_form, const std::vector&lt; SAMRAI::solv::RobinBcCoefStrategy&lt; NDIM &gt; * &gt; &amp;bc_coefs)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~INSCollocatedPPMConvectiveOperator</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_n_s_collocated_p_p_m_convective_operator.html</anchorfile>
      <anchor>a0f942dd514cf6e0c3d877efcfa58410b</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>applyConvectiveOperator</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_n_s_collocated_p_p_m_convective_operator.html</anchorfile>
      <anchor>a747e5f8248c23494d069fd5d9515454a</anchor>
      <arglist>(int U_idx, int N_idx)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>initializeOperatorState</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_n_s_collocated_p_p_m_convective_operator.html</anchorfile>
      <anchor>abf531ef32eda58427e5a74daffb43f9c</anchor>
      <arglist>(const SAMRAI::solv::SAMRAIVectorReal&lt; NDIM, double &gt; &amp;in, const SAMRAI::solv::SAMRAIVectorReal&lt; NDIM, double &gt; &amp;out)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>deallocateOperatorState</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_n_s_collocated_p_p_m_convective_operator.html</anchorfile>
      <anchor>a3f78f1a071591d8e3ed67047e47a1635</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" static="yes">
      <type>static SAMRAI::tbox::Pointer&lt; ConvectiveOperator &gt;</type>
      <name>allocate_operator</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_n_s_collocated_p_p_m_convective_operator.html</anchorfile>
      <anchor>a50841bc5c9d4823a814ad6cb7b600235</anchor>
      <arglist>(const std::string &amp;object_name, SAMRAI::tbox::Pointer&lt; SAMRAI::tbox::Database &gt; input_db, ConvectiveDifferencingType difference_form, const std::vector&lt; SAMRAI::solv::RobinBcCoefStrategy&lt; NDIM &gt; * &gt; &amp;bc_coefs)</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>IBAMR::INSCollocatedVelocityBcCoef</name>
    <filename>class_i_b_a_m_r_1_1_i_n_s_collocated_velocity_bc_coef.html</filename>
    <base>IBAMR::StokesBcCoefStrategy</base>
    <member kind="function">
      <type></type>
      <name>INSCollocatedVelocityBcCoef</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_n_s_collocated_velocity_bc_coef.html</anchorfile>
      <anchor>a641639b14a964c853dfa7cae31a7364f</anchor>
      <arglist>(unsigned int comp_idx, const INSCollocatedHierarchyIntegrator *fluid_solver, const std::vector&lt; SAMRAI::solv::RobinBcCoefStrategy&lt; NDIM &gt; * &gt; &amp;bc_coefs, TractionBcType traction_bc_type, bool homogeneous_bc=false)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~INSCollocatedVelocityBcCoef</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_n_s_collocated_velocity_bc_coef.html</anchorfile>
      <anchor>a12b3c0e8d552f6defba619c93342cf2c</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setPhysicalBcCoefs</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_n_s_collocated_velocity_bc_coef.html</anchorfile>
      <anchor>a104081be162ab5dca80f9821ccf8bc88</anchor>
      <arglist>(const std::vector&lt; SAMRAI::solv::RobinBcCoefStrategy&lt; NDIM &gt; * &gt; &amp;bc_coefs)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setSolutionTime</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_n_s_collocated_velocity_bc_coef.html</anchorfile>
      <anchor>a897952e84739206038510983b8205c05</anchor>
      <arglist>(double solution_time)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setTimeInterval</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_n_s_collocated_velocity_bc_coef.html</anchorfile>
      <anchor>ae65ca4308d80d40c6b05be4b14a79e7f</anchor>
      <arglist>(double current_time, double new_time)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setStokesSpecifications</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_n_s_collocated_velocity_bc_coef.html</anchorfile>
      <anchor>aa50c3b7f694acf739ec079960709a592</anchor>
      <arglist>(const StokesSpecifications *problem_coefs)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setTargetVelocityPatchDataIndex</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_n_s_collocated_velocity_bc_coef.html</anchorfile>
      <anchor>a3d27d2880cb25c0b3478cfff5e49379f</anchor>
      <arglist>(int u_target_data_idx)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>clearTargetVelocityPatchDataIndex</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_n_s_collocated_velocity_bc_coef.html</anchorfile>
      <anchor>a9c7e864ba50993ca24a141ba90575a9b</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setTargetPressurePatchDataIndex</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_n_s_collocated_velocity_bc_coef.html</anchorfile>
      <anchor>a72be05a207888d2d528ae1d08f244c01</anchor>
      <arglist>(int p_target_data_idx)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>clearTargetPressurePatchDataIndex</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_n_s_collocated_velocity_bc_coef.html</anchorfile>
      <anchor>a6eadc68cf5123c611c251167d29c8899</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setTargetPatchDataIndex</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_n_s_collocated_velocity_bc_coef.html</anchorfile>
      <anchor>a6c5752c88f015a4470d3767371ef3b07</anchor>
      <arglist>(int target_idx)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>clearTargetPatchDataIndex</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_n_s_collocated_velocity_bc_coef.html</anchorfile>
      <anchor>a6561252d0fae696fc02489990fa23e9f</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setHomogeneousBc</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_n_s_collocated_velocity_bc_coef.html</anchorfile>
      <anchor>a6aed87bf07f67dd19771bc4b13fb0a4e</anchor>
      <arglist>(bool homogeneous_bc)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setBcCoefs</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_n_s_collocated_velocity_bc_coef.html</anchorfile>
      <anchor>ad2380e6de7083a4daf7e80402a70bcf6</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::ArrayData&lt; NDIM, double &gt; &gt; &amp;acoef_data, SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::ArrayData&lt; NDIM, double &gt; &gt; &amp;bcoef_data, SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::ArrayData&lt; NDIM, double &gt; &gt; &amp;gcoef_data, const SAMRAI::tbox::Pointer&lt; SAMRAI::hier::Variable&lt; NDIM &gt; &gt; &amp;variable, const SAMRAI::hier::Patch&lt; NDIM &gt; &amp;patch, const SAMRAI::hier::BoundaryBox&lt; NDIM &gt; &amp;bdry_box, double fill_time=0.0) const </arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>StokesBcCoefStrategy</name>
      <anchorfile>class_i_b_a_m_r_1_1_stokes_bc_coef_strategy.html</anchorfile>
      <anchor>aaff128acb8e64aa0360af6c55fa2f3f8</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~StokesBcCoefStrategy</name>
      <anchorfile>class_i_b_a_m_r_1_1_stokes_bc_coef_strategy.html</anchorfile>
      <anchor>a77428a4b35a43d9cd5aecf961b770735</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>setTractionBcType</name>
      <anchorfile>class_i_b_a_m_r_1_1_stokes_bc_coef_strategy.html</anchorfile>
      <anchor>a752b9043b0e49b5029b3ef321de573d4</anchor>
      <arglist>(TractionBcType bc_type)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual TractionBcType</type>
      <name>getTractionBcType</name>
      <anchorfile>class_i_b_a_m_r_1_1_stokes_bc_coef_strategy.html</anchorfile>
      <anchor>a2368851d57d5e4877986b33d03d8b532</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>int</type>
      <name>d_u_target_data_idx</name>
      <anchorfile>class_i_b_a_m_r_1_1_stokes_bc_coef_strategy.html</anchorfile>
      <anchor>a1350b69c6805490c214fa248cb6756e3</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>IBAMR::INSHierarchyIntegrator</name>
    <filename>class_i_b_a_m_r_1_1_i_n_s_hierarchy_integrator.html</filename>
    <base>IBTK::HierarchyIntegrator</base>
    <member kind="function" virtualness="pure">
      <type>virtual SAMRAI::tbox::Pointer&lt; ConvectiveOperator &gt;</type>
      <name>getConvectiveOperator</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_n_s_hierarchy_integrator.html</anchorfile>
      <anchor>a703790dbc2ad134a576f641d23b1480b</anchor>
      <arglist>()=0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual SAMRAI::tbox::Pointer&lt; IBTK::PoissonSolver &gt;</type>
      <name>getVelocitySubdomainSolver</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_n_s_hierarchy_integrator.html</anchorfile>
      <anchor>abf2a8cc2af8318078f76455e6e800f9b</anchor>
      <arglist>()=0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual SAMRAI::tbox::Pointer&lt; IBTK::PoissonSolver &gt;</type>
      <name>getPressureSubdomainSolver</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_n_s_hierarchy_integrator.html</anchorfile>
      <anchor>aed1f2bc8dc283798087e1aefe291c359</anchor>
      <arglist>()=0</arglist>
    </member>
    <member kind="function" protection="protected" virtualness="pure">
      <type>virtual double</type>
      <name>getStableTimestep</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_n_s_hierarchy_integrator.html</anchorfile>
      <anchor>afc0ac0b0a89903257b422b8755a72829</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::hier::Patch&lt; NDIM &gt; &gt; patch) const =0</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>IBAMR::INSIntermediateVelocityBcCoef</name>
    <filename>class_i_b_a_m_r_1_1_i_n_s_intermediate_velocity_bc_coef.html</filename>
    <base>IBTK::ExtendedRobinBcCoefStrategy</base>
    <member kind="function">
      <type></type>
      <name>INSIntermediateVelocityBcCoef</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_n_s_intermediate_velocity_bc_coef.html</anchorfile>
      <anchor>a8b24d121f655acf5761f83b0d89c2575</anchor>
      <arglist>(int comp_idx, const std::vector&lt; SAMRAI::solv::RobinBcCoefStrategy&lt; NDIM &gt; * &gt; &amp;bc_coefs, bool homogeneous_bc=false)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~INSIntermediateVelocityBcCoef</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_n_s_intermediate_velocity_bc_coef.html</anchorfile>
      <anchor>aea44f1ae7f3e37c9961ec2453d9df37c</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setPhysicalBcCoefs</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_n_s_intermediate_velocity_bc_coef.html</anchorfile>
      <anchor>a14d68457881972f4f289f42947331daa</anchor>
      <arglist>(const std::vector&lt; SAMRAI::solv::RobinBcCoefStrategy&lt; NDIM &gt; * &gt; &amp;bc_coefs)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setSolutionTime</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_n_s_intermediate_velocity_bc_coef.html</anchorfile>
      <anchor>a71be0f4958f335ff81c3c87b797a364a</anchor>
      <arglist>(double solution_time)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setTimeInterval</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_n_s_intermediate_velocity_bc_coef.html</anchorfile>
      <anchor>ab843c7c7327f7c4089eae7564dd3e55d</anchor>
      <arglist>(double current_time, double new_time)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setTargetPatchDataIndex</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_n_s_intermediate_velocity_bc_coef.html</anchorfile>
      <anchor>a0a4f7d9f4afa8b7202ea92ea88e34390</anchor>
      <arglist>(int target_idx)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>clearTargetPatchDataIndex</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_n_s_intermediate_velocity_bc_coef.html</anchorfile>
      <anchor>ad257fc089bce4600b7b55ce8ec464050</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setHomogeneousBc</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_n_s_intermediate_velocity_bc_coef.html</anchorfile>
      <anchor>abfc7063c836a3fce642973c84b44474b</anchor>
      <arglist>(bool homogeneous_bc)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setBcCoefs</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_n_s_intermediate_velocity_bc_coef.html</anchorfile>
      <anchor>a39829f27afb364e613014a311395193b</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::ArrayData&lt; NDIM, double &gt; &gt; &amp;acoef_data, SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::ArrayData&lt; NDIM, double &gt; &gt; &amp;bcoef_data, SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::ArrayData&lt; NDIM, double &gt; &gt; &amp;gcoef_data, const SAMRAI::tbox::Pointer&lt; SAMRAI::hier::Variable&lt; NDIM &gt; &gt; &amp;variable, const SAMRAI::hier::Patch&lt; NDIM &gt; &amp;patch, const SAMRAI::hier::BoundaryBox&lt; NDIM &gt; &amp;bdry_box, double fill_time=0.0) const </arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>IBAMR::INSProjectionBcCoef</name>
    <filename>class_i_b_a_m_r_1_1_i_n_s_projection_bc_coef.html</filename>
    <base>IBTK::ExtendedRobinBcCoefStrategy</base>
    <member kind="function">
      <type></type>
      <name>INSProjectionBcCoef</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_n_s_projection_bc_coef.html</anchorfile>
      <anchor>a4c9feffa1be139b7aa30a92f6b876787</anchor>
      <arglist>(const std::vector&lt; SAMRAI::solv::RobinBcCoefStrategy&lt; NDIM &gt; * &gt; &amp;bc_coefs, bool homogeneous_bc=false)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~INSProjectionBcCoef</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_n_s_projection_bc_coef.html</anchorfile>
      <anchor>af9f50e358906ace73ff8cfea83850b3d</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setPhysicalBcCoefs</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_n_s_projection_bc_coef.html</anchorfile>
      <anchor>a8f1288a11c358d4ed707d9eb5f383da4</anchor>
      <arglist>(const std::vector&lt; SAMRAI::solv::RobinBcCoefStrategy&lt; NDIM &gt; * &gt; &amp;bc_coefs)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setSolutionTime</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_n_s_projection_bc_coef.html</anchorfile>
      <anchor>a619af15f3f7810ca2c1685b2ed0e24de</anchor>
      <arglist>(double solution_time)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setTimeInterval</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_n_s_projection_bc_coef.html</anchorfile>
      <anchor>ab8dcfbc9d430045abaef44c29dfc442b</anchor>
      <arglist>(double current_time, double new_time)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setTargetPatchDataIndex</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_n_s_projection_bc_coef.html</anchorfile>
      <anchor>a1456fe045347e5b29e044f9b4a4c66a0</anchor>
      <arglist>(int target_idx)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>clearTargetPatchDataIndex</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_n_s_projection_bc_coef.html</anchorfile>
      <anchor>a55b5356d5d279a0454272296cd57d952</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setHomogeneousBc</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_n_s_projection_bc_coef.html</anchorfile>
      <anchor>a4769fa68985e544e083889b431f767c1</anchor>
      <arglist>(bool homogeneous_bc)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setBcCoefs</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_n_s_projection_bc_coef.html</anchorfile>
      <anchor>a3a2cbfcca31b16bc83c171bb2f1b9cc3</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::ArrayData&lt; NDIM, double &gt; &gt; &amp;acoef_data, SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::ArrayData&lt; NDIM, double &gt; &gt; &amp;bcoef_data, SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::ArrayData&lt; NDIM, double &gt; &gt; &amp;gcoef_data, const SAMRAI::tbox::Pointer&lt; SAMRAI::hier::Variable&lt; NDIM &gt; &gt; &amp;variable, const SAMRAI::hier::Patch&lt; NDIM &gt; &amp;patch, const SAMRAI::hier::BoundaryBox&lt; NDIM &gt; &amp;bdry_box, double fill_time=0.0) const </arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>IBAMR::INSStaggeredCenteredConvectiveOperator</name>
    <filename>class_i_b_a_m_r_1_1_i_n_s_staggered_centered_convective_operator.html</filename>
    <base>IBAMR::ConvectiveOperator</base>
    <member kind="function">
      <type></type>
      <name>INSStaggeredCenteredConvectiveOperator</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_n_s_staggered_centered_convective_operator.html</anchorfile>
      <anchor>a0bc49da56cdc111fb0623d33b21740c1</anchor>
      <arglist>(const std::string &amp;object_name, SAMRAI::tbox::Pointer&lt; SAMRAI::tbox::Database &gt; input_db, ConvectiveDifferencingType difference_form, const std::vector&lt; SAMRAI::solv::RobinBcCoefStrategy&lt; NDIM &gt; * &gt; &amp;bc_coefs)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~INSStaggeredCenteredConvectiveOperator</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_n_s_staggered_centered_convective_operator.html</anchorfile>
      <anchor>acf12b919e9058af74c527e0226a23ea4</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>applyConvectiveOperator</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_n_s_staggered_centered_convective_operator.html</anchorfile>
      <anchor>ad3069a3d934106e5f51374801cdd8351</anchor>
      <arglist>(int U_idx, int N_idx)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>initializeOperatorState</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_n_s_staggered_centered_convective_operator.html</anchorfile>
      <anchor>a7e1a6186997f3307182b53aa1e933b31</anchor>
      <arglist>(const SAMRAI::solv::SAMRAIVectorReal&lt; NDIM, double &gt; &amp;in, const SAMRAI::solv::SAMRAIVectorReal&lt; NDIM, double &gt; &amp;out)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>deallocateOperatorState</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_n_s_staggered_centered_convective_operator.html</anchorfile>
      <anchor>abc765301e250542a1f55eac47d6530a6</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" static="yes">
      <type>static SAMRAI::tbox::Pointer&lt; ConvectiveOperator &gt;</type>
      <name>allocate_operator</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_n_s_staggered_centered_convective_operator.html</anchorfile>
      <anchor>ae13d8cb7f07aae24ea9d57e6dad290fb</anchor>
      <arglist>(const std::string &amp;object_name, SAMRAI::tbox::Pointer&lt; SAMRAI::tbox::Database &gt; input_db, ConvectiveDifferencingType difference_form, const std::vector&lt; SAMRAI::solv::RobinBcCoefStrategy&lt; NDIM &gt; * &gt; &amp;bc_coefs)</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>IBAMR::INSStaggeredConvectiveOperatorManager</name>
    <filename>class_i_b_a_m_r_1_1_i_n_s_staggered_convective_operator_manager.html</filename>
    <member kind="typedef">
      <type>SAMRAI::tbox::Pointer&lt; ConvectiveOperator &gt;(*</type>
      <name>OperatorMaker</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_n_s_staggered_convective_operator_manager.html</anchorfile>
      <anchor>a7dee1181a9087bfd878f2550f82934ed</anchor>
      <arglist>)(const std::string &amp;operator_object_name, SAMRAI::tbox::Pointer&lt; SAMRAI::tbox::Database &gt; input_db, ConvectiveDifferencingType difference_form, const std::vector&lt; SAMRAI::solv::RobinBcCoefStrategy&lt; NDIM &gt; * &gt; &amp;bc_coefs)</arglist>
    </member>
    <member kind="function">
      <type>SAMRAI::tbox::Pointer&lt; ConvectiveOperator &gt;</type>
      <name>allocateOperator</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_n_s_staggered_convective_operator_manager.html</anchorfile>
      <anchor>a11d43b5aa708d7ed8811b943b6d86c90</anchor>
      <arglist>(const std::string &amp;operator_type, const std::string &amp;operator_object_name, SAMRAI::tbox::Pointer&lt; SAMRAI::tbox::Database &gt; input_db, ConvectiveDifferencingType difference_form, const std::vector&lt; SAMRAI::solv::RobinBcCoefStrategy&lt; NDIM &gt; * &gt; &amp;bc_coefs) const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>registerOperatorFactoryFunction</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_n_s_staggered_convective_operator_manager.html</anchorfile>
      <anchor>a11bd2dc94e8b9ddf05e28d28dbe26e39</anchor>
      <arglist>(const std::string &amp;operator_type, OperatorMaker operator_maker)</arglist>
    </member>
    <member kind="function" static="yes">
      <type>static INSStaggeredConvectiveOperatorManager *</type>
      <name>getManager</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_n_s_staggered_convective_operator_manager.html</anchorfile>
      <anchor>a0e6a6d965d0ea70af4e74437247f58fc</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" static="yes">
      <type>static void</type>
      <name>freeManager</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_n_s_staggered_convective_operator_manager.html</anchorfile>
      <anchor>a96da0fadea0175e4fd134fad19699347</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="variable" static="yes">
      <type>static const std::string</type>
      <name>DEFAULT</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_n_s_staggered_convective_operator_manager.html</anchorfile>
      <anchor>af850806f7a8111a512fd6c2dabd4ccef</anchor>
      <arglist></arglist>
    </member>
    <member kind="function" protection="protected">
      <type></type>
      <name>INSStaggeredConvectiveOperatorManager</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_n_s_staggered_convective_operator_manager.html</anchorfile>
      <anchor>a3bff7387caae1571615e9b54cf62169a</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" protection="protected">
      <type></type>
      <name>~INSStaggeredConvectiveOperatorManager</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_n_s_staggered_convective_operator_manager.html</anchorfile>
      <anchor>af742f9007af1e27c32169a7c28fc0932</anchor>
      <arglist>()</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>IBAMR::INSStaggeredHierarchyIntegrator</name>
    <filename>class_i_b_a_m_r_1_1_i_n_s_staggered_hierarchy_integrator.html</filename>
    <base>IBAMR::INSHierarchyIntegrator</base>
    <member kind="function">
      <type></type>
      <name>INSStaggeredHierarchyIntegrator</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_n_s_staggered_hierarchy_integrator.html</anchorfile>
      <anchor>acf8546864f2ea219c15f63189bdf2469</anchor>
      <arglist>(const std::string &amp;object_name, SAMRAI::tbox::Pointer&lt; SAMRAI::tbox::Database &gt; input_db, bool register_for_restart=true)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~INSStaggeredHierarchyIntegrator</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_n_s_staggered_hierarchy_integrator.html</anchorfile>
      <anchor>aba7094b4e7caa45fe2bc48dddce522b0</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>SAMRAI::tbox::Pointer&lt; ConvectiveOperator &gt;</type>
      <name>getConvectiveOperator</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_n_s_staggered_hierarchy_integrator.html</anchorfile>
      <anchor>ad38bc8c169820ba7d6a4c980d6d0fa5f</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>SAMRAI::tbox::Pointer&lt; IBTK::PoissonSolver &gt;</type>
      <name>getVelocitySubdomainSolver</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_n_s_staggered_hierarchy_integrator.html</anchorfile>
      <anchor>a7b06a9c51c3da1ff30d4c05fcefc7ac5</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>SAMRAI::tbox::Pointer&lt; IBTK::PoissonSolver &gt;</type>
      <name>getPressureSubdomainSolver</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_n_s_staggered_hierarchy_integrator.html</anchorfile>
      <anchor>a47e38cb09ebceab15beea31e607c19ab</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setStokesSolver</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_n_s_staggered_hierarchy_integrator.html</anchorfile>
      <anchor>a5cb40f28d581a523d687177e595c3cfc</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; StaggeredStokesSolver &gt; stokes_solver)</arglist>
    </member>
    <member kind="function">
      <type>SAMRAI::tbox::Pointer&lt; StaggeredStokesSolver &gt;</type>
      <name>getStokesSolver</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_n_s_staggered_hierarchy_integrator.html</anchorfile>
      <anchor>a96bc4e5e24949cfe617dcbdbfe1cf659</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>initializeHierarchyIntegrator</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_n_s_staggered_hierarchy_integrator.html</anchorfile>
      <anchor>ab4bde9547da273f75b6afd1cb5c2ef01</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::hier::PatchHierarchy&lt; NDIM &gt; &gt; hierarchy, SAMRAI::tbox::Pointer&lt; SAMRAI::mesh::GriddingAlgorithm&lt; NDIM &gt; &gt; gridding_alg)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>initializePatchHierarchy</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_n_s_staggered_hierarchy_integrator.html</anchorfile>
      <anchor>ad865e17c4dbcc87ec0412b4328291e77</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::hier::PatchHierarchy&lt; NDIM &gt; &gt; hierarchy, SAMRAI::tbox::Pointer&lt; SAMRAI::mesh::GriddingAlgorithm&lt; NDIM &gt; &gt; gridding_alg)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>preprocessIntegrateHierarchy</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_n_s_staggered_hierarchy_integrator.html</anchorfile>
      <anchor>a41accb4f51887c674aa500c3e9b67aca</anchor>
      <arglist>(double current_time, double new_time, int num_cycles=1)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>integrateHierarchy</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_n_s_staggered_hierarchy_integrator.html</anchorfile>
      <anchor>aaad59d7f901cafe3d880178b257fb31d</anchor>
      <arglist>(double current_time, double new_time, int cycle_num=0)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>postprocessIntegrateHierarchy</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_n_s_staggered_hierarchy_integrator.html</anchorfile>
      <anchor>adcc1f9aa986162d0622ecacbb1910df4</anchor>
      <arglist>(double current_time, double new_time, bool skip_synchronize_new_state_data, int num_cycles=1)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>regridHierarchy</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_n_s_staggered_hierarchy_integrator.html</anchorfile>
      <anchor>ad3f1b8a845d8edd72de007efe8ac918a</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setupSolverVectors</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_n_s_staggered_hierarchy_integrator.html</anchorfile>
      <anchor>a26f70e125e127adf0c9284b8cce0d949</anchor>
      <arglist>(const SAMRAI::tbox::Pointer&lt; SAMRAI::solv::SAMRAIVectorReal&lt; NDIM, double &gt; &gt; &amp;sol_vec, const SAMRAI::tbox::Pointer&lt; SAMRAI::solv::SAMRAIVectorReal&lt; NDIM, double &gt; &gt; &amp;rhs_vec, double current_time, double new_time, int cycle_num)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>resetSolverVectors</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_n_s_staggered_hierarchy_integrator.html</anchorfile>
      <anchor>a83d3facec56cccf756e6fba19435ca55</anchor>
      <arglist>(const SAMRAI::tbox::Pointer&lt; SAMRAI::solv::SAMRAIVectorReal&lt; NDIM, double &gt; &gt; &amp;sol_vec, const SAMRAI::tbox::Pointer&lt; SAMRAI::solv::SAMRAIVectorReal&lt; NDIM, double &gt; &gt; &amp;rhs_vec, double current_time, double new_time, int cycle_num)</arglist>
    </member>
    <member kind="function" protection="protected">
      <type>double</type>
      <name>getStableTimestep</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_n_s_staggered_hierarchy_integrator.html</anchorfile>
      <anchor>a1fed2aa98c6f0b84cccbb3519d4b1e86</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::hier::Patch&lt; NDIM &gt; &gt; patch) const </arglist>
    </member>
    <member kind="function" protection="protected">
      <type>void</type>
      <name>initializeLevelDataSpecialized</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_n_s_staggered_hierarchy_integrator.html</anchorfile>
      <anchor>a01573418a3188ad653d2c120319d1c20</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::hier::BasePatchHierarchy&lt; NDIM &gt; &gt; hierarchy, int level_number, double init_data_time, bool can_be_refined, bool initial_time, SAMRAI::tbox::Pointer&lt; SAMRAI::hier::BasePatchLevel&lt; NDIM &gt; &gt; old_level, bool allocate_data)</arglist>
    </member>
    <member kind="function" protection="protected">
      <type>void</type>
      <name>resetHierarchyConfigurationSpecialized</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_n_s_staggered_hierarchy_integrator.html</anchorfile>
      <anchor>ac3b20a049111e144bac05b9a09486b3e</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::hier::BasePatchHierarchy&lt; NDIM &gt; &gt; hierarchy, int coarsest_level, int finest_level)</arglist>
    </member>
    <member kind="function" protection="protected">
      <type>void</type>
      <name>applyGradientDetectorSpecialized</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_n_s_staggered_hierarchy_integrator.html</anchorfile>
      <anchor>a76be9d741080fb6e31a76bc4108a8b96</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::hier::BasePatchHierarchy&lt; NDIM &gt; &gt; hierarchy, int level_number, double error_data_time, int tag_index, bool initial_time, bool uses_richardson_extrapolation_too)</arglist>
    </member>
    <member kind="function" protection="protected">
      <type>void</type>
      <name>setupPlotDataSpecialized</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_n_s_staggered_hierarchy_integrator.html</anchorfile>
      <anchor>a8a59dbc34dc948ccce0c5d78776797c2</anchor>
      <arglist>()</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>IBAMR::INSStaggeredPPMConvectiveOperator</name>
    <filename>class_i_b_a_m_r_1_1_i_n_s_staggered_p_p_m_convective_operator.html</filename>
    <base>IBAMR::ConvectiveOperator</base>
    <member kind="function">
      <type></type>
      <name>INSStaggeredPPMConvectiveOperator</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_n_s_staggered_p_p_m_convective_operator.html</anchorfile>
      <anchor>a856b4b581a7ed1d148bbc6480b22f369</anchor>
      <arglist>(const std::string &amp;object_name, SAMRAI::tbox::Pointer&lt; SAMRAI::tbox::Database &gt; input_db, ConvectiveDifferencingType difference_form, const std::vector&lt; SAMRAI::solv::RobinBcCoefStrategy&lt; NDIM &gt; * &gt; &amp;bc_coefs)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~INSStaggeredPPMConvectiveOperator</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_n_s_staggered_p_p_m_convective_operator.html</anchorfile>
      <anchor>aa2efaf195a2a501d060102c941c44bae</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>applyConvectiveOperator</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_n_s_staggered_p_p_m_convective_operator.html</anchorfile>
      <anchor>a652cdc3281aeb413d689c1ce8a4bf6ff</anchor>
      <arglist>(int U_idx, int N_idx)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>initializeOperatorState</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_n_s_staggered_p_p_m_convective_operator.html</anchorfile>
      <anchor>a492d6ccd695007cb7f52d829632f11da</anchor>
      <arglist>(const SAMRAI::solv::SAMRAIVectorReal&lt; NDIM, double &gt; &amp;in, const SAMRAI::solv::SAMRAIVectorReal&lt; NDIM, double &gt; &amp;out)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>deallocateOperatorState</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_n_s_staggered_p_p_m_convective_operator.html</anchorfile>
      <anchor>ab3bb7158762a5daaec55281c6f837c4c</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" static="yes">
      <type>static SAMRAI::tbox::Pointer&lt; ConvectiveOperator &gt;</type>
      <name>allocate_operator</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_n_s_staggered_p_p_m_convective_operator.html</anchorfile>
      <anchor>a7fe241ae831cec3815815169fab3843a</anchor>
      <arglist>(const std::string &amp;object_name, SAMRAI::tbox::Pointer&lt; SAMRAI::tbox::Database &gt; input_db, ConvectiveDifferencingType difference_form, const std::vector&lt; SAMRAI::solv::RobinBcCoefStrategy&lt; NDIM &gt; * &gt; &amp;bc_coefs)</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>IBAMR::INSStaggeredPressureBcCoef</name>
    <filename>class_i_b_a_m_r_1_1_i_n_s_staggered_pressure_bc_coef.html</filename>
    <base>IBAMR::StokesBcCoefStrategy</base>
    <member kind="function">
      <type></type>
      <name>INSStaggeredPressureBcCoef</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_n_s_staggered_pressure_bc_coef.html</anchorfile>
      <anchor>abfe0a4376e5cf9a4cf8ec76e2e4022e4</anchor>
      <arglist>(const INSStaggeredHierarchyIntegrator *fluid_solver, const std::vector&lt; SAMRAI::solv::RobinBcCoefStrategy&lt; NDIM &gt; * &gt; &amp;bc_coefs, TractionBcType traction_bc_type, bool homogeneous_bc=false)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~INSStaggeredPressureBcCoef</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_n_s_staggered_pressure_bc_coef.html</anchorfile>
      <anchor>a353572fa90ed0bd7cc9cf4593bf7a877</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setPhysicalBcCoefs</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_n_s_staggered_pressure_bc_coef.html</anchorfile>
      <anchor>a917bba774c932382709f2eda9823063c</anchor>
      <arglist>(const std::vector&lt; SAMRAI::solv::RobinBcCoefStrategy&lt; NDIM &gt; * &gt; &amp;bc_coefs)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setSolutionTime</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_n_s_staggered_pressure_bc_coef.html</anchorfile>
      <anchor>ab6e24a2cf274d2f970e27b083467b844</anchor>
      <arglist>(double solution_time)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setTimeInterval</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_n_s_staggered_pressure_bc_coef.html</anchorfile>
      <anchor>ae3821ab3668ef6643b22c265689af0ea</anchor>
      <arglist>(double current_time, double new_time)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setStokesSpecifications</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_n_s_staggered_pressure_bc_coef.html</anchorfile>
      <anchor>a35384393fe28663a7485010e15e5d26e</anchor>
      <arglist>(const StokesSpecifications *problem_coefs)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setTargetVelocityPatchDataIndex</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_n_s_staggered_pressure_bc_coef.html</anchorfile>
      <anchor>adb7be51f2150208cf73ebacc71c2a26b</anchor>
      <arglist>(int u_target_data_idx)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>clearTargetVelocityPatchDataIndex</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_n_s_staggered_pressure_bc_coef.html</anchorfile>
      <anchor>a38473db577a3c5f3176a1cdbdc418dc3</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setTargetPressurePatchDataIndex</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_n_s_staggered_pressure_bc_coef.html</anchorfile>
      <anchor>ad1dabf001feb45c41193f649d63ec7a8</anchor>
      <arglist>(int p_target_data_idx)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>clearTargetPressurePatchDataIndex</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_n_s_staggered_pressure_bc_coef.html</anchorfile>
      <anchor>aa95e27f45a8a8513534b77bbce3c5914</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setTargetPatchDataIndex</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_n_s_staggered_pressure_bc_coef.html</anchorfile>
      <anchor>ab631806e8108d4a82f0091e323825f90</anchor>
      <arglist>(int target_idx)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>clearTargetPatchDataIndex</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_n_s_staggered_pressure_bc_coef.html</anchorfile>
      <anchor>aaa1f1c7dbdf8ab786245e6ac9fabf6bc</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setHomogeneousBc</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_n_s_staggered_pressure_bc_coef.html</anchorfile>
      <anchor>ae84e2a75a0555ea939b7b110f619b78d</anchor>
      <arglist>(bool homogeneous_bc)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setBcCoefs</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_n_s_staggered_pressure_bc_coef.html</anchorfile>
      <anchor>a9a5aaf01ab071662e64b1165b5d8d057</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::ArrayData&lt; NDIM, double &gt; &gt; &amp;acoef_data, SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::ArrayData&lt; NDIM, double &gt; &gt; &amp;bcoef_data, SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::ArrayData&lt; NDIM, double &gt; &gt; &amp;gcoef_data, const SAMRAI::tbox::Pointer&lt; SAMRAI::hier::Variable&lt; NDIM &gt; &gt; &amp;variable, const SAMRAI::hier::Patch&lt; NDIM &gt; &amp;patch, const SAMRAI::hier::BoundaryBox&lt; NDIM &gt; &amp;bdry_box, double fill_time=0.0) const </arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>IBAMR::INSStaggeredStabilizedPPMConvectiveOperator</name>
    <filename>class_i_b_a_m_r_1_1_i_n_s_staggered_stabilized_p_p_m_convective_operator.html</filename>
    <base>IBAMR::ConvectiveOperator</base>
    <member kind="function">
      <type></type>
      <name>INSStaggeredStabilizedPPMConvectiveOperator</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_n_s_staggered_stabilized_p_p_m_convective_operator.html</anchorfile>
      <anchor>a8ee99b99dec09a0229db25f176b9d855</anchor>
      <arglist>(const std::string &amp;object_name, SAMRAI::tbox::Pointer&lt; SAMRAI::tbox::Database &gt; input_db, ConvectiveDifferencingType difference_form, const std::vector&lt; SAMRAI::solv::RobinBcCoefStrategy&lt; NDIM &gt; * &gt; &amp;bc_coefs)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~INSStaggeredStabilizedPPMConvectiveOperator</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_n_s_staggered_stabilized_p_p_m_convective_operator.html</anchorfile>
      <anchor>adafdc6a80bd5bc701f654cc0d5b5f2b9</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>applyConvectiveOperator</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_n_s_staggered_stabilized_p_p_m_convective_operator.html</anchorfile>
      <anchor>abe3b5ed7e096076c050a440a8d5e9dca</anchor>
      <arglist>(int U_idx, int N_idx)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>initializeOperatorState</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_n_s_staggered_stabilized_p_p_m_convective_operator.html</anchorfile>
      <anchor>a03dc469ca35155dbf29baa05c162897d</anchor>
      <arglist>(const SAMRAI::solv::SAMRAIVectorReal&lt; NDIM, double &gt; &amp;in, const SAMRAI::solv::SAMRAIVectorReal&lt; NDIM, double &gt; &amp;out)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>deallocateOperatorState</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_n_s_staggered_stabilized_p_p_m_convective_operator.html</anchorfile>
      <anchor>a80a29ca37366decdad33c329edba5a3d</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" static="yes">
      <type>static SAMRAI::tbox::Pointer&lt; ConvectiveOperator &gt;</type>
      <name>allocate_operator</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_n_s_staggered_stabilized_p_p_m_convective_operator.html</anchorfile>
      <anchor>ab9de0ff17ed35f5dba0c46f93447b6a1</anchor>
      <arglist>(const std::string &amp;object_name, SAMRAI::tbox::Pointer&lt; SAMRAI::tbox::Database &gt; input_db, ConvectiveDifferencingType difference_form, const std::vector&lt; SAMRAI::solv::RobinBcCoefStrategy&lt; NDIM &gt; * &gt; &amp;bc_coefs)</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>IBAMR::INSStaggeredStochasticForcing</name>
    <filename>class_i_b_a_m_r_1_1_i_n_s_staggered_stochastic_forcing.html</filename>
    <base>IBTK::CartGridFunction</base>
    <member kind="function">
      <type></type>
      <name>INSStaggeredStochasticForcing</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_n_s_staggered_stochastic_forcing.html</anchorfile>
      <anchor>ae2c77bc80fc9f31956a44746152bd908</anchor>
      <arglist>(const std::string &amp;object_name, SAMRAI::tbox::Pointer&lt; SAMRAI::tbox::Database &gt; input_db, const INSStaggeredHierarchyIntegrator *fluid_solver)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~INSStaggeredStochasticForcing</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_n_s_staggered_stochastic_forcing.html</anchorfile>
      <anchor>a488e699f6409c36fe820abd1928a3559</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>isTimeDependent</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_n_s_staggered_stochastic_forcing.html</anchorfile>
      <anchor>aa52cf2f7be3b19cd331f39b39eb8d591</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setDataOnPatchHierarchy</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_n_s_staggered_stochastic_forcing.html</anchorfile>
      <anchor>a1177e777a61ad00406d6a95396a3cae8</anchor>
      <arglist>(const int data_idx, SAMRAI::tbox::Pointer&lt; SAMRAI::hier::Variable&lt; NDIM &gt; &gt; var, SAMRAI::tbox::Pointer&lt; SAMRAI::hier::PatchHierarchy&lt; NDIM &gt; &gt; hierarchy, const double data_time, const bool initial_time=false, const int coarsest_ln=-1, const int finest_ln=-1)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setDataOnPatch</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_n_s_staggered_stochastic_forcing.html</anchorfile>
      <anchor>ac44bf6e0e1f6b31c0817c0a364b3f54d</anchor>
      <arglist>(const int data_idx, SAMRAI::tbox::Pointer&lt; SAMRAI::hier::Variable&lt; NDIM &gt; &gt; var, SAMRAI::tbox::Pointer&lt; SAMRAI::hier::Patch&lt; NDIM &gt; &gt; patch, const double data_time, const bool initial_time=false, SAMRAI::tbox::Pointer&lt; SAMRAI::hier::PatchLevel&lt; NDIM &gt; &gt; patch_level=SAMRAI::tbox::Pointer&lt; SAMRAI::hier::PatchLevel&lt; NDIM &gt; &gt;(NULL))</arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>std::string</type>
      <name>d_object_name</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_n_s_staggered_stochastic_forcing.html</anchorfile>
      <anchor>a674fa78856c9a8d3a36099a03e3d2b3e</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>IBAMR::INSStaggeredUpwindConvectiveOperator</name>
    <filename>class_i_b_a_m_r_1_1_i_n_s_staggered_upwind_convective_operator.html</filename>
    <base>IBAMR::ConvectiveOperator</base>
    <member kind="function">
      <type></type>
      <name>INSStaggeredUpwindConvectiveOperator</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_n_s_staggered_upwind_convective_operator.html</anchorfile>
      <anchor>ab8cb0afcd82b936221e389e662f7e462</anchor>
      <arglist>(const std::string &amp;object_name, SAMRAI::tbox::Pointer&lt; SAMRAI::tbox::Database &gt; input_db, ConvectiveDifferencingType difference_form, const std::vector&lt; SAMRAI::solv::RobinBcCoefStrategy&lt; NDIM &gt; * &gt; &amp;bc_coefs)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~INSStaggeredUpwindConvectiveOperator</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_n_s_staggered_upwind_convective_operator.html</anchorfile>
      <anchor>a68a98856c66ea6f38eb39b4fd021e835</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>applyConvectiveOperator</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_n_s_staggered_upwind_convective_operator.html</anchorfile>
      <anchor>ae2562c122c693eb451572db16b58688a</anchor>
      <arglist>(int U_idx, int N_idx)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>initializeOperatorState</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_n_s_staggered_upwind_convective_operator.html</anchorfile>
      <anchor>aa6ae9d406cec5fe69a4efa644836ce3a</anchor>
      <arglist>(const SAMRAI::solv::SAMRAIVectorReal&lt; NDIM, double &gt; &amp;in, const SAMRAI::solv::SAMRAIVectorReal&lt; NDIM, double &gt; &amp;out)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>deallocateOperatorState</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_n_s_staggered_upwind_convective_operator.html</anchorfile>
      <anchor>a4634ca2c6d576433025ca88a2ea9a868</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" static="yes">
      <type>static SAMRAI::tbox::Pointer&lt; ConvectiveOperator &gt;</type>
      <name>allocate_operator</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_n_s_staggered_upwind_convective_operator.html</anchorfile>
      <anchor>aceabd7a50c33d4a7f6a4783cf2a34fa2</anchor>
      <arglist>(const std::string &amp;object_name, SAMRAI::tbox::Pointer&lt; SAMRAI::tbox::Database &gt; input_db, ConvectiveDifferencingType difference_form, const std::vector&lt; SAMRAI::solv::RobinBcCoefStrategy&lt; NDIM &gt; * &gt; &amp;bc_coefs)</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>IBAMR::INSStaggeredVelocityBcCoef</name>
    <filename>class_i_b_a_m_r_1_1_i_n_s_staggered_velocity_bc_coef.html</filename>
    <base>IBAMR::StokesBcCoefStrategy</base>
    <member kind="function">
      <type></type>
      <name>INSStaggeredVelocityBcCoef</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_n_s_staggered_velocity_bc_coef.html</anchorfile>
      <anchor>a44701e6953e70486ff61f991ee91867c</anchor>
      <arglist>(unsigned int comp_idx, const INSStaggeredHierarchyIntegrator *fluid_solver, const std::vector&lt; SAMRAI::solv::RobinBcCoefStrategy&lt; NDIM &gt; * &gt; &amp;bc_coefs, TractionBcType traction_bc_type, bool homogeneous_bc=false)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~INSStaggeredVelocityBcCoef</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_n_s_staggered_velocity_bc_coef.html</anchorfile>
      <anchor>a9ee12617cd02b5633706c3c813b547bb</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setPhysicalBcCoefs</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_n_s_staggered_velocity_bc_coef.html</anchorfile>
      <anchor>a7b85eb166bb187ea14306adb314a1134</anchor>
      <arglist>(const std::vector&lt; SAMRAI::solv::RobinBcCoefStrategy&lt; NDIM &gt; * &gt; &amp;bc_coefs)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setSolutionTime</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_n_s_staggered_velocity_bc_coef.html</anchorfile>
      <anchor>a4a87596c8ae9f26397ca7a38b73ddaee</anchor>
      <arglist>(double solution_time)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setTimeInterval</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_n_s_staggered_velocity_bc_coef.html</anchorfile>
      <anchor>a2f1c8f819daa3d6b5e776b4dc07b1b8c</anchor>
      <arglist>(double current_time, double new_time)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setStokesSpecifications</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_n_s_staggered_velocity_bc_coef.html</anchorfile>
      <anchor>a5040a3aa8716ecaf2940fe9f95231841</anchor>
      <arglist>(const StokesSpecifications *problem_coefs)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setTargetVelocityPatchDataIndex</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_n_s_staggered_velocity_bc_coef.html</anchorfile>
      <anchor>ad4e60075c65aad59d628cfffc6c01e05</anchor>
      <arglist>(int u_target_data_idx)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>clearTargetVelocityPatchDataIndex</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_n_s_staggered_velocity_bc_coef.html</anchorfile>
      <anchor>a06f7a1502b81954bd76e3ffc0cc136ec</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setTargetPressurePatchDataIndex</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_n_s_staggered_velocity_bc_coef.html</anchorfile>
      <anchor>a0287c65d4ace6999d553bc9af72f60a7</anchor>
      <arglist>(int p_target_data_idx)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>clearTargetPressurePatchDataIndex</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_n_s_staggered_velocity_bc_coef.html</anchorfile>
      <anchor>ae66c794caa2ab1524ec3beccf91b7590</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setTargetPatchDataIndex</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_n_s_staggered_velocity_bc_coef.html</anchorfile>
      <anchor>ad61187cfdc8b1410515264b5e9ce51a5</anchor>
      <arglist>(int target_idx)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>clearTargetPatchDataIndex</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_n_s_staggered_velocity_bc_coef.html</anchorfile>
      <anchor>a7d1f64f9791f870f02ef26494bc561ce</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setHomogeneousBc</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_n_s_staggered_velocity_bc_coef.html</anchorfile>
      <anchor>a8b13aa336ebc2ace38019fd30be60bf3</anchor>
      <arglist>(bool homogeneous_bc)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setBcCoefs</name>
      <anchorfile>class_i_b_a_m_r_1_1_i_n_s_staggered_velocity_bc_coef.html</anchorfile>
      <anchor>ac0640f4790045e05631115e1ac32cf54</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::ArrayData&lt; NDIM, double &gt; &gt; &amp;acoef_data, SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::ArrayData&lt; NDIM, double &gt; &gt; &amp;bcoef_data, SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::ArrayData&lt; NDIM, double &gt; &gt; &amp;gcoef_data, const SAMRAI::tbox::Pointer&lt; SAMRAI::hier::Variable&lt; NDIM &gt; &gt; &amp;variable, const SAMRAI::hier::Patch&lt; NDIM &gt; &amp;patch, const SAMRAI::hier::BoundaryBox&lt; NDIM &gt; &amp;bdry_box, double fill_time=0.0) const </arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>IBAMR::KrylovLinearSolverStaggeredStokesSolverInterface</name>
    <filename>class_i_b_a_m_r_1_1_krylov_linear_solver_staggered_stokes_solver_interface.html</filename>
    <base>IBAMR::StaggeredStokesSolver</base>
    <member kind="function">
      <type></type>
      <name>KrylovLinearSolverStaggeredStokesSolverInterface</name>
      <anchorfile>class_i_b_a_m_r_1_1_krylov_linear_solver_staggered_stokes_solver_interface.html</anchorfile>
      <anchor>a14cc855d50fa1e3d963b59bb96fe70ea</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~KrylovLinearSolverStaggeredStokesSolverInterface</name>
      <anchorfile>class_i_b_a_m_r_1_1_krylov_linear_solver_staggered_stokes_solver_interface.html</anchorfile>
      <anchor>afa31b27f10938e3be143614f2ec32cb2</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setVelocityPoissonSpecifications</name>
      <anchorfile>class_i_b_a_m_r_1_1_krylov_linear_solver_staggered_stokes_solver_interface.html</anchorfile>
      <anchor>a802625cdf4b929dd9f1fe7fef14e7b0a</anchor>
      <arglist>(const SAMRAI::solv::PoissonSpecifications &amp;U_problem_coefs)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setPhysicalBcCoefs</name>
      <anchorfile>class_i_b_a_m_r_1_1_krylov_linear_solver_staggered_stokes_solver_interface.html</anchorfile>
      <anchor>af458338612cc5b840ac47e87e9601194</anchor>
      <arglist>(const std::vector&lt; SAMRAI::solv::RobinBcCoefStrategy&lt; NDIM &gt; * &gt; &amp;U_bc_coefs, SAMRAI::solv::RobinBcCoefStrategy&lt; NDIM &gt; *P_bc_coef)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setPhysicalBoundaryHelper</name>
      <anchorfile>class_i_b_a_m_r_1_1_krylov_linear_solver_staggered_stokes_solver_interface.html</anchorfile>
      <anchor>a8ff6ef4b91ab42ed099474edc703007f</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; StaggeredStokesPhysicalBoundaryHelper &gt; bc_helper)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>StaggeredStokesSolver</name>
      <anchorfile>class_i_b_a_m_r_1_1_staggered_stokes_solver.html</anchorfile>
      <anchor>aa669dcc798a052a8554f4bf46fc886d6</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~StaggeredStokesSolver</name>
      <anchorfile>class_i_b_a_m_r_1_1_staggered_stokes_solver.html</anchorfile>
      <anchor>ae58acab154fa31a966517bea5ac8e719</anchor>
      <arglist>()</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>IBAMR::MaterialPointSpec</name>
    <filename>class_i_b_a_m_r_1_1_material_point_spec.html</filename>
    <base>IBTK::Streamable</base>
    <member kind="function">
      <type></type>
      <name>MaterialPointSpec</name>
      <anchorfile>class_i_b_a_m_r_1_1_material_point_spec.html</anchorfile>
      <anchor>a5a68d382fe1242812cbb51002930978b</anchor>
      <arglist>(int point_idx=-1, double weight=0.0, libMesh::subdomain_id_type subdomain_id=0, const std::vector&lt; double &gt; &amp;internal_vars=std::vector&lt; double &gt;())</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~MaterialPointSpec</name>
      <anchorfile>class_i_b_a_m_r_1_1_material_point_spec.html</anchorfile>
      <anchor>a221535e61eec85649f93ad913aeeac3c</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>const int &amp;</type>
      <name>getPointIndex</name>
      <anchorfile>class_i_b_a_m_r_1_1_material_point_spec.html</anchorfile>
      <anchor>aa68866d73442ca2682e46f231cd2ca3d</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>int &amp;</type>
      <name>getPointIndex</name>
      <anchorfile>class_i_b_a_m_r_1_1_material_point_spec.html</anchorfile>
      <anchor>ad1a97e54e4000f87f2e8c7208028104a</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>const double &amp;</type>
      <name>getWeight</name>
      <anchorfile>class_i_b_a_m_r_1_1_material_point_spec.html</anchorfile>
      <anchor>ac241ea1e12fcdc5b013c1d0a7be6dc38</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>double &amp;</type>
      <name>getWeight</name>
      <anchorfile>class_i_b_a_m_r_1_1_material_point_spec.html</anchorfile>
      <anchor>a5a840225665d97560a80fda4c4f7965d</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>const libMesh::subdomain_id_type &amp;</type>
      <name>getSubdomainId</name>
      <anchorfile>class_i_b_a_m_r_1_1_material_point_spec.html</anchorfile>
      <anchor>a8f295c3d6fbf7ab5464c27b9a62190a3</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>libMesh::subdomain_id_type &amp;</type>
      <name>getSubdomainId</name>
      <anchorfile>class_i_b_a_m_r_1_1_material_point_spec.html</anchorfile>
      <anchor>a576452930a2ccffd09d541626cac3784</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>const std::vector&lt; double &gt; &amp;</type>
      <name>getInternalVariables</name>
      <anchorfile>class_i_b_a_m_r_1_1_material_point_spec.html</anchorfile>
      <anchor>a8a06a748beadb481c5b02091bc682dc4</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>std::vector&lt; double &gt; &amp;</type>
      <name>getInternalVariables</name>
      <anchorfile>class_i_b_a_m_r_1_1_material_point_spec.html</anchorfile>
      <anchor>aabb405bc4342e2db22c704a93f1e1b2f</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>getStreamableClassID</name>
      <anchorfile>class_i_b_a_m_r_1_1_material_point_spec.html</anchorfile>
      <anchor>aaaf6ae7ae25a17a8d5b457ed4e29b9d1</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>size_t</type>
      <name>getDataStreamSize</name>
      <anchorfile>class_i_b_a_m_r_1_1_material_point_spec.html</anchorfile>
      <anchor>a9cdef0addd6811a17bc33b63fd974c24</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>packStream</name>
      <anchorfile>class_i_b_a_m_r_1_1_material_point_spec.html</anchorfile>
      <anchor>a5c42b2e27d4e574a526e14a60c7bf8cf</anchor>
      <arglist>(SAMRAI::tbox::AbstractStream &amp;stream)</arglist>
    </member>
    <member kind="function" static="yes">
      <type>static void</type>
      <name>registerWithStreamableManager</name>
      <anchorfile>class_i_b_a_m_r_1_1_material_point_spec.html</anchorfile>
      <anchor>abb03abdf9e902229e34c7d8af3309aae</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" static="yes">
      <type>static bool</type>
      <name>getIsRegisteredWithStreamableManager</name>
      <anchorfile>class_i_b_a_m_r_1_1_material_point_spec.html</anchorfile>
      <anchor>a1aba6b0541b03ea45b7ff983ebd3eb48</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="variable" static="yes">
      <type>static int</type>
      <name>STREAMABLE_CLASS_ID</name>
      <anchorfile>class_i_b_a_m_r_1_1_material_point_spec.html</anchorfile>
      <anchor>a4b5ab970968c1e185e6f0b391eac4c65</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>IBAMR::PenaltyIBMethod</name>
    <filename>class_i_b_a_m_r_1_1_penalty_i_b_method.html</filename>
    <base>IBAMR::IBMethod</base>
    <member kind="function">
      <type></type>
      <name>PenaltyIBMethod</name>
      <anchorfile>class_i_b_a_m_r_1_1_penalty_i_b_method.html</anchorfile>
      <anchor>ad1a00a9b0f4ba9a42d72e4cee7bc505d</anchor>
      <arglist>(const std::string &amp;object_name, SAMRAI::tbox::Pointer&lt; SAMRAI::tbox::Database &gt; input_db, bool register_for_restart=true)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~PenaltyIBMethod</name>
      <anchorfile>class_i_b_a_m_r_1_1_penalty_i_b_method.html</anchorfile>
      <anchor>a5c0b019c2b0a94bea96473bc2e14f988</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>preprocessIntegrateData</name>
      <anchorfile>class_i_b_a_m_r_1_1_penalty_i_b_method.html</anchorfile>
      <anchor>afe2ab4dba45f2beaa87c6af9eb43df77</anchor>
      <arglist>(double current_time, double new_time, int num_cycles)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>postprocessIntegrateData</name>
      <anchorfile>class_i_b_a_m_r_1_1_penalty_i_b_method.html</anchorfile>
      <anchor>ade4f6f67597c75a879e65b3d04e1b6ed</anchor>
      <arglist>(double current_time, double new_time, int num_cycles)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>eulerStep</name>
      <anchorfile>class_i_b_a_m_r_1_1_penalty_i_b_method.html</anchorfile>
      <anchor>afb85f3286c55d0a4227f47bb3d419927</anchor>
      <arglist>(double current_time, double new_time)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>midpointStep</name>
      <anchorfile>class_i_b_a_m_r_1_1_penalty_i_b_method.html</anchorfile>
      <anchor>a15f7e01f940b6a53a9270e5e13eab97e</anchor>
      <arglist>(double current_time, double new_time)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>trapezoidalStep</name>
      <anchorfile>class_i_b_a_m_r_1_1_penalty_i_b_method.html</anchorfile>
      <anchor>abe42f890d78e4f3d2639fb04d024cb52</anchor>
      <arglist>(double current_time, double new_time)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>computeLagrangianForce</name>
      <anchorfile>class_i_b_a_m_r_1_1_penalty_i_b_method.html</anchorfile>
      <anchor>a7f54eea00feb988bd8e7b8dc00ce9623</anchor>
      <arglist>(double data_time)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>initializePatchHierarchy</name>
      <anchorfile>class_i_b_a_m_r_1_1_penalty_i_b_method.html</anchorfile>
      <anchor>abdbd6218e1149eb1d0571d801c1a612f</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::hier::PatchHierarchy&lt; NDIM &gt; &gt; hierarchy, SAMRAI::tbox::Pointer&lt; SAMRAI::mesh::GriddingAlgorithm&lt; NDIM &gt; &gt; gridding_alg, int u_data_idx, const std::vector&lt; SAMRAI::tbox::Pointer&lt; SAMRAI::xfer::CoarsenSchedule&lt; NDIM &gt; &gt; &gt; &amp;u_synch_scheds, const std::vector&lt; SAMRAI::tbox::Pointer&lt; SAMRAI::xfer::RefineSchedule&lt; NDIM &gt; &gt; &gt; &amp;u_ghost_fill_scheds, int integrator_step, double init_data_time, bool initial_time)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>putToDatabase</name>
      <anchorfile>class_i_b_a_m_r_1_1_penalty_i_b_method.html</anchorfile>
      <anchor>a728be9a015501026e462e46e1ec4a440</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::tbox::Database &gt; db)</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>IBAMR::PETScKrylovStaggeredStokesSolver</name>
    <filename>class_i_b_a_m_r_1_1_p_e_t_sc_krylov_staggered_stokes_solver.html</filename>
    <base>IBTK::PETScKrylovLinearSolver</base>
    <base>IBAMR::KrylovLinearSolverStaggeredStokesSolverInterface</base>
    <member kind="function">
      <type></type>
      <name>PETScKrylovStaggeredStokesSolver</name>
      <anchorfile>class_i_b_a_m_r_1_1_p_e_t_sc_krylov_staggered_stokes_solver.html</anchorfile>
      <anchor>acb8f668d99781926fb204fb72782dfba</anchor>
      <arglist>(const std::string &amp;object_name, SAMRAI::tbox::Pointer&lt; SAMRAI::tbox::Database &gt; input_db, const std::string &amp;default_options_prefix)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~PETScKrylovStaggeredStokesSolver</name>
      <anchorfile>class_i_b_a_m_r_1_1_p_e_t_sc_krylov_staggered_stokes_solver.html</anchorfile>
      <anchor>a3f4a4c3abb2b11da8904d04dab0be999</anchor>
      <arglist>()</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>IBAMR::RNG</name>
    <filename>class_i_b_a_m_r_1_1_r_n_g.html</filename>
  </compound>
  <compound kind="class">
    <name>IBAMR::SpongeLayerForceFunction</name>
    <filename>class_i_b_a_m_r_1_1_sponge_layer_force_function.html</filename>
    <base>IBTK::CartGridFunction</base>
    <member kind="function">
      <type></type>
      <name>SpongeLayerForceFunction</name>
      <anchorfile>class_i_b_a_m_r_1_1_sponge_layer_force_function.html</anchorfile>
      <anchor>a8eceb521df5a299b025ebe9f37e9a2bb</anchor>
      <arglist>(const std::string &amp;object_name, SAMRAI::tbox::Pointer&lt; SAMRAI::tbox::Database &gt; input_db, const INSHierarchyIntegrator *fluid_solver, SAMRAI::tbox::Pointer&lt; SAMRAI::geom::CartesianGridGeometry&lt; NDIM &gt; &gt; grid_geometry)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~SpongeLayerForceFunction</name>
      <anchorfile>class_i_b_a_m_r_1_1_sponge_layer_force_function.html</anchorfile>
      <anchor>ab53c531fd2f80ea889e332820fdb34f3</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>isTimeDependent</name>
      <anchorfile>class_i_b_a_m_r_1_1_sponge_layer_force_function.html</anchorfile>
      <anchor>ac8436da55dd0f79222c2c88ce04f6437</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setDataOnPatch</name>
      <anchorfile>class_i_b_a_m_r_1_1_sponge_layer_force_function.html</anchorfile>
      <anchor>a5913d83e2fe031d7a5033aeaf0fed2fb</anchor>
      <arglist>(int data_idx, SAMRAI::tbox::Pointer&lt; SAMRAI::hier::Variable&lt; NDIM &gt; &gt; var, SAMRAI::tbox::Pointer&lt; SAMRAI::hier::Patch&lt; NDIM &gt; &gt; patch, double data_time, bool initial_time=false, SAMRAI::tbox::Pointer&lt; SAMRAI::hier::PatchLevel&lt; NDIM &gt; &gt; level=SAMRAI::tbox::Pointer&lt; SAMRAI::hier::PatchLevel&lt; NDIM &gt; &gt;(NULL))</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>IBAMR::StaggeredStokesBlockFactorizationPreconditioner</name>
    <filename>class_i_b_a_m_r_1_1_staggered_stokes_block_factorization_preconditioner.html</filename>
    <base>IBAMR::StaggeredStokesBlockPreconditioner</base>
    <member kind="function">
      <type></type>
      <name>StaggeredStokesBlockFactorizationPreconditioner</name>
      <anchorfile>class_i_b_a_m_r_1_1_staggered_stokes_block_factorization_preconditioner.html</anchorfile>
      <anchor>a5dc1732688bf79a0f3fa306c3c99fe3a</anchor>
      <arglist>(const std::string &amp;object_name, SAMRAI::tbox::Pointer&lt; SAMRAI::tbox::Database &gt; input_db, const std::string &amp;default_options_prefix)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~StaggeredStokesBlockFactorizationPreconditioner</name>
      <anchorfile>class_i_b_a_m_r_1_1_staggered_stokes_block_factorization_preconditioner.html</anchorfile>
      <anchor>a045c064e29ce3403a5f542e28bd0ba90</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>solveSystem</name>
      <anchorfile>class_i_b_a_m_r_1_1_staggered_stokes_block_factorization_preconditioner.html</anchorfile>
      <anchor>a06b2dbbdb5c23bae696ce009f390d5cf</anchor>
      <arglist>(SAMRAI::solv::SAMRAIVectorReal&lt; NDIM, double &gt; &amp;x, SAMRAI::solv::SAMRAIVectorReal&lt; NDIM, double &gt; &amp;b)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>initializeSolverState</name>
      <anchorfile>class_i_b_a_m_r_1_1_staggered_stokes_block_factorization_preconditioner.html</anchorfile>
      <anchor>ad8c97021776935e09bd5bb31aaf41107</anchor>
      <arglist>(const SAMRAI::solv::SAMRAIVectorReal&lt; NDIM, double &gt; &amp;x, const SAMRAI::solv::SAMRAIVectorReal&lt; NDIM, double &gt; &amp;b)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>deallocateSolverState</name>
      <anchorfile>class_i_b_a_m_r_1_1_staggered_stokes_block_factorization_preconditioner.html</anchorfile>
      <anchor>a6f99ec28739da2c42e2c8ffa22dfbf04</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setInitialGuessNonzero</name>
      <anchorfile>class_i_b_a_m_r_1_1_staggered_stokes_block_factorization_preconditioner.html</anchorfile>
      <anchor>a4a986550e0627091a270d9d9956a6e5f</anchor>
      <arglist>(bool initial_guess_nonzero=true)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setMaxIterations</name>
      <anchorfile>class_i_b_a_m_r_1_1_staggered_stokes_block_factorization_preconditioner.html</anchorfile>
      <anchor>ac73c51dd9b9ed80fd085a40d6efd6995</anchor>
      <arglist>(int max_iterations)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>StaggeredStokesBlockPreconditioner</name>
      <anchorfile>class_i_b_a_m_r_1_1_staggered_stokes_block_preconditioner.html</anchorfile>
      <anchor>afd3bbc58f09cd74dd483f534b15dfe9e</anchor>
      <arglist>(bool needs_velocity_solver, bool needs_pressure_solver)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~StaggeredStokesBlockPreconditioner</name>
      <anchorfile>class_i_b_a_m_r_1_1_staggered_stokes_block_preconditioner.html</anchorfile>
      <anchor>a8dabc679546532b7ea8369b367cf49a4</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual bool</type>
      <name>needsVelocitySubdomainSolver</name>
      <anchorfile>class_i_b_a_m_r_1_1_staggered_stokes_block_preconditioner.html</anchorfile>
      <anchor>a036ed995877cb79f1197beb0bfe6b1c3</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>setVelocitySubdomainSolver</name>
      <anchorfile>class_i_b_a_m_r_1_1_staggered_stokes_block_preconditioner.html</anchorfile>
      <anchor>a81c3158ac62fcd9e4a6602bbb56132f5</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; IBTK::PoissonSolver &gt; velocity_solver)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setVelocityPoissonSpecifications</name>
      <anchorfile>class_i_b_a_m_r_1_1_staggered_stokes_block_preconditioner.html</anchorfile>
      <anchor>a659350d6a2290b77f3f4c6d8d8d9265f</anchor>
      <arglist>(const SAMRAI::solv::PoissonSpecifications &amp;U_problem_coefs)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual bool</type>
      <name>needsPressureSubdomainSolver</name>
      <anchorfile>class_i_b_a_m_r_1_1_staggered_stokes_block_preconditioner.html</anchorfile>
      <anchor>aa3fa32b7f519ac6a4703b686c698c52a</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>setPressureSubdomainSolver</name>
      <anchorfile>class_i_b_a_m_r_1_1_staggered_stokes_block_preconditioner.html</anchorfile>
      <anchor>a646f15a73d72be51c3e72fe0425c0cb4</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; IBTK::PoissonSolver &gt; pressure_solver)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>setPressurePoissonSpecifications</name>
      <anchorfile>class_i_b_a_m_r_1_1_staggered_stokes_block_preconditioner.html</anchorfile>
      <anchor>ade767ec1008b14e45c985d78a88fd918</anchor>
      <arglist>(const SAMRAI::solv::PoissonSpecifications &amp;P_problem_coefs)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>setPhysicalBcCoefs</name>
      <anchorfile>class_i_b_a_m_r_1_1_staggered_stokes_block_preconditioner.html</anchorfile>
      <anchor>a6b0817b1ce133f90800521913d4f19cc</anchor>
      <arglist>(const std::vector&lt; SAMRAI::solv::RobinBcCoefStrategy&lt; NDIM &gt; * &gt; &amp;U_bc_coefs, SAMRAI::solv::RobinBcCoefStrategy&lt; NDIM &gt; *P_bc_coef)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>initializeSolverState</name>
      <anchorfile>class_i_b_a_m_r_1_1_staggered_stokes_block_preconditioner.html</anchorfile>
      <anchor>ab4363051afe51d8ea8efcf92b52558c6</anchor>
      <arglist>(const SAMRAI::solv::SAMRAIVectorReal&lt; NDIM, double &gt; &amp;x, const SAMRAI::solv::SAMRAIVectorReal&lt; NDIM, double &gt; &amp;b)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>deallocateSolverState</name>
      <anchorfile>class_i_b_a_m_r_1_1_staggered_stokes_block_preconditioner.html</anchorfile>
      <anchor>a81f71c7a982b6ef5fc2563b752d0f02f</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>setPhysicalBoundaryHelper</name>
      <anchorfile>class_i_b_a_m_r_1_1_staggered_stokes_solver.html</anchorfile>
      <anchor>ac8637f24e6d3fb4f148aa836a50017e6</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; StaggeredStokesPhysicalBoundaryHelper &gt; bc_helper)</arglist>
    </member>
    <member kind="function" static="yes">
      <type>static SAMRAI::tbox::Pointer&lt; StaggeredStokesSolver &gt;</type>
      <name>allocate_solver</name>
      <anchorfile>class_i_b_a_m_r_1_1_staggered_stokes_block_factorization_preconditioner.html</anchorfile>
      <anchor>ae5209ff60bd87d60126c97dd2ae5a02d</anchor>
      <arglist>(const std::string &amp;object_name, SAMRAI::tbox::Pointer&lt; SAMRAI::tbox::Database &gt; input_db, const std::string &amp;default_options_prefix)</arglist>
    </member>
    <member kind="function" protection="protected">
      <type>void</type>
      <name>correctNullspace</name>
      <anchorfile>class_i_b_a_m_r_1_1_staggered_stokes_block_preconditioner.html</anchorfile>
      <anchor>aefc4ba72d457068dc5e2c17536fe2798</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::solv::SAMRAIVectorReal&lt; NDIM, double &gt; &gt; U_vec, SAMRAI::tbox::Pointer&lt; SAMRAI::solv::SAMRAIVectorReal&lt; NDIM, double &gt; &gt; P_vec)</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>IBAMR::StaggeredStokesBlockPreconditioner</name>
    <filename>class_i_b_a_m_r_1_1_staggered_stokes_block_preconditioner.html</filename>
    <base>IBTK::LinearSolver</base>
    <base>IBAMR::StaggeredStokesSolver</base>
  </compound>
  <compound kind="class">
    <name>IBAMR::StaggeredStokesBoxRelaxationFACOperator</name>
    <filename>class_i_b_a_m_r_1_1_staggered_stokes_box_relaxation_f_a_c_operator.html</filename>
    <base>IBAMR::StaggeredStokesFACPreconditionerStrategy</base>
    <member kind="function">
      <type></type>
      <name>StaggeredStokesBoxRelaxationFACOperator</name>
      <anchorfile>class_i_b_a_m_r_1_1_staggered_stokes_box_relaxation_f_a_c_operator.html</anchorfile>
      <anchor>a3fd8efb850506a3a88a0668976da7010</anchor>
      <arglist>(const std::string &amp;object_name, SAMRAI::tbox::Pointer&lt; SAMRAI::tbox::Database &gt; input_db, const std::string &amp;default_options_prefix)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~StaggeredStokesBoxRelaxationFACOperator</name>
      <anchorfile>class_i_b_a_m_r_1_1_staggered_stokes_box_relaxation_f_a_c_operator.html</anchorfile>
      <anchor>a7dfdd5aa74b7c5fb9a8aa6fbbb32be53</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>smoothError</name>
      <anchorfile>class_i_b_a_m_r_1_1_staggered_stokes_box_relaxation_f_a_c_operator.html</anchorfile>
      <anchor>a2d975294aed47187fa692717bf832d69</anchor>
      <arglist>(SAMRAI::solv::SAMRAIVectorReal&lt; NDIM, double &gt; &amp;error, const SAMRAI::solv::SAMRAIVectorReal&lt; NDIM, double &gt; &amp;residual, int level_num, int num_sweeps, bool performing_pre_sweeps, bool performing_post_sweeps)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>StaggeredStokesFACPreconditionerStrategy</name>
      <anchorfile>class_i_b_a_m_r_1_1_staggered_stokes_f_a_c_preconditioner_strategy.html</anchorfile>
      <anchor>a529acb25bb6df1d59c32fa893bf47c09</anchor>
      <arglist>(const std::string &amp;object_name, int ghost_cell_width, SAMRAI::tbox::Pointer&lt; SAMRAI::tbox::Database &gt; input_db, const std::string &amp;default_options_prefix)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~StaggeredStokesFACPreconditionerStrategy</name>
      <anchorfile>class_i_b_a_m_r_1_1_staggered_stokes_f_a_c_preconditioner_strategy.html</anchorfile>
      <anchor>a0aeb91b1ca85f7dba643eec60d16ac1f</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>setVelocityPoissonSpecifications</name>
      <anchorfile>class_i_b_a_m_r_1_1_staggered_stokes_f_a_c_preconditioner_strategy.html</anchorfile>
      <anchor>a0dfc5f828dd17aa8371d135d505122e5</anchor>
      <arglist>(const SAMRAI::solv::PoissonSpecifications &amp;U_problem_coefs)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>setPhysicalBcCoefs</name>
      <anchorfile>class_i_b_a_m_r_1_1_staggered_stokes_f_a_c_preconditioner_strategy.html</anchorfile>
      <anchor>a534d1d9047d2060a9269eb4dde4b1474</anchor>
      <arglist>(const std::vector&lt; SAMRAI::solv::RobinBcCoefStrategy&lt; NDIM &gt; * &gt; &amp;U_bc_coefs, SAMRAI::solv::RobinBcCoefStrategy&lt; NDIM &gt; *P_bc_coef)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>setPhysicalBoundaryHelper</name>
      <anchorfile>class_i_b_a_m_r_1_1_staggered_stokes_f_a_c_preconditioner_strategy.html</anchorfile>
      <anchor>a39ea03feaf000ba3366a5b846ef54ae3</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; StaggeredStokesPhysicalBoundaryHelper &gt; bc_helper)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setResetLevels</name>
      <anchorfile>class_i_b_a_m_r_1_1_staggered_stokes_f_a_c_preconditioner_strategy.html</anchorfile>
      <anchor>ad6fef9a860b8c730f16f884e620f10d7</anchor>
      <arglist>(int coarsest_ln, int finest_ln)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setSmootherType</name>
      <anchorfile>class_i_b_a_m_r_1_1_staggered_stokes_f_a_c_preconditioner_strategy.html</anchorfile>
      <anchor>a775abb92030c7c87395ba3ab4fde67c9</anchor>
      <arglist>(const std::string &amp;smoother_type)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setCoarseSolverType</name>
      <anchorfile>class_i_b_a_m_r_1_1_staggered_stokes_f_a_c_preconditioner_strategy.html</anchorfile>
      <anchor>ae82f8885b880c16778e66ee726070fb1</anchor>
      <arglist>(const std::string &amp;coarse_solver_type)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setCoarseSolverMaxIterations</name>
      <anchorfile>class_i_b_a_m_r_1_1_staggered_stokes_f_a_c_preconditioner_strategy.html</anchorfile>
      <anchor>ab3773227be34e7ae516b86edc2b7eb8c</anchor>
      <arglist>(int coarse_solver_max_iterations)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setCoarseSolverAbsoluteTolerance</name>
      <anchorfile>class_i_b_a_m_r_1_1_staggered_stokes_f_a_c_preconditioner_strategy.html</anchorfile>
      <anchor>a21b08f3f06e4449bc8ba1a8d6d429eb6</anchor>
      <arglist>(double coarse_solver_abs_residual_tol)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setCoarseSolverRelativeTolerance</name>
      <anchorfile>class_i_b_a_m_r_1_1_staggered_stokes_f_a_c_preconditioner_strategy.html</anchorfile>
      <anchor>a776b117d38412a138600f22fc8cd71b7</anchor>
      <arglist>(double coarse_solver_rel_residual_tol)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setProlongationMethods</name>
      <anchorfile>class_i_b_a_m_r_1_1_staggered_stokes_f_a_c_preconditioner_strategy.html</anchorfile>
      <anchor>aba97f31f103c67a35d08b76719c15663</anchor>
      <arglist>(const std::string &amp;U_prolongation_method, const std::string &amp;P_prolongation_method)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setRestrictionMethods</name>
      <anchorfile>class_i_b_a_m_r_1_1_staggered_stokes_f_a_c_preconditioner_strategy.html</anchorfile>
      <anchor>a411abaf3d5fc095474a4bf08c166529f</anchor>
      <arglist>(const std::string &amp;U_restriction_method, const std::string &amp;P_restriction_method)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>restrictResidual</name>
      <anchorfile>class_i_b_a_m_r_1_1_staggered_stokes_f_a_c_preconditioner_strategy.html</anchorfile>
      <anchor>a3b964942e0f1bb90aa8cb046ced203a4</anchor>
      <arglist>(const SAMRAI::solv::SAMRAIVectorReal&lt; NDIM, double &gt; &amp;src, SAMRAI::solv::SAMRAIVectorReal&lt; NDIM, double &gt; &amp;dst, int dst_ln)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>prolongError</name>
      <anchorfile>class_i_b_a_m_r_1_1_staggered_stokes_f_a_c_preconditioner_strategy.html</anchorfile>
      <anchor>adb536365b3655c6f2662490201a08d2b</anchor>
      <arglist>(const SAMRAI::solv::SAMRAIVectorReal&lt; NDIM, double &gt; &amp;src, SAMRAI::solv::SAMRAIVectorReal&lt; NDIM, double &gt; &amp;dst, int dst_ln)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>prolongErrorAndCorrect</name>
      <anchorfile>class_i_b_a_m_r_1_1_staggered_stokes_f_a_c_preconditioner_strategy.html</anchorfile>
      <anchor>a3984848e99b98a5e04b747b7dfa4f141</anchor>
      <arglist>(const SAMRAI::solv::SAMRAIVectorReal&lt; NDIM, double &gt; &amp;src, SAMRAI::solv::SAMRAIVectorReal&lt; NDIM, double &gt; &amp;dst, int dst_ln)</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>solveCoarsestLevel</name>
      <anchorfile>class_i_b_a_m_r_1_1_staggered_stokes_f_a_c_preconditioner_strategy.html</anchorfile>
      <anchor>adf678458a4516da676da263c591b4220</anchor>
      <arglist>(SAMRAI::solv::SAMRAIVectorReal&lt; NDIM, double &gt; &amp;error, const SAMRAI::solv::SAMRAIVectorReal&lt; NDIM, double &gt; &amp;residual, int coarsest_ln)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>computeResidual</name>
      <anchorfile>class_i_b_a_m_r_1_1_staggered_stokes_f_a_c_preconditioner_strategy.html</anchorfile>
      <anchor>aed24e75e261e2cbe9cbe1c3a71f13e87</anchor>
      <arglist>(SAMRAI::solv::SAMRAIVectorReal&lt; NDIM, double &gt; &amp;residual, const SAMRAI::solv::SAMRAIVectorReal&lt; NDIM, double &gt; &amp;solution, const SAMRAI::solv::SAMRAIVectorReal&lt; NDIM, double &gt; &amp;rhs, int coarsest_level_num, int finest_level_num)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>initializeOperatorState</name>
      <anchorfile>class_i_b_a_m_r_1_1_staggered_stokes_f_a_c_preconditioner_strategy.html</anchorfile>
      <anchor>add53c0d4eff28ec905e92503be72ce8e</anchor>
      <arglist>(const SAMRAI::solv::SAMRAIVectorReal&lt; NDIM, double &gt; &amp;solution, const SAMRAI::solv::SAMRAIVectorReal&lt; NDIM, double &gt; &amp;rhs)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>deallocateOperatorState</name>
      <anchorfile>class_i_b_a_m_r_1_1_staggered_stokes_f_a_c_preconditioner_strategy.html</anchorfile>
      <anchor>a978d9439481ca0cd839293b534e730fd</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>allocateScratchData</name>
      <anchorfile>class_i_b_a_m_r_1_1_staggered_stokes_f_a_c_preconditioner_strategy.html</anchorfile>
      <anchor>aa9676806a07802bc6ef3c23172d0a5de</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>deallocateScratchData</name>
      <anchorfile>class_i_b_a_m_r_1_1_staggered_stokes_f_a_c_preconditioner_strategy.html</anchorfile>
      <anchor>a5ccb176a56217024170ff9c648019c1c</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" protection="protected">
      <type>void</type>
      <name>initializeOperatorStateSpecialized</name>
      <anchorfile>class_i_b_a_m_r_1_1_staggered_stokes_box_relaxation_f_a_c_operator.html</anchorfile>
      <anchor>af4308a694402e1fe7d4869c4b14546f4</anchor>
      <arglist>(const SAMRAI::solv::SAMRAIVectorReal&lt; NDIM, double &gt; &amp;solution, const SAMRAI::solv::SAMRAIVectorReal&lt; NDIM, double &gt; &amp;rhs, int coarsest_reset_ln, int finest_reset_ln)</arglist>
    </member>
    <member kind="function" protection="protected">
      <type>void</type>
      <name>deallocateOperatorStateSpecialized</name>
      <anchorfile>class_i_b_a_m_r_1_1_staggered_stokes_box_relaxation_f_a_c_operator.html</anchorfile>
      <anchor>a5d26f1ea5bd853fee4e9a87ef46aa9e3</anchor>
      <arglist>(int coarsest_reset_ln, int finest_reset_ln)</arglist>
    </member>
    <member kind="function" protection="protected">
      <type>void</type>
      <name>xeqScheduleProlongation</name>
      <anchorfile>class_i_b_a_m_r_1_1_staggered_stokes_f_a_c_preconditioner_strategy.html</anchorfile>
      <anchor>a0463f734ae67eaa8a1af1074f183fdae</anchor>
      <arglist>(const std::pair&lt; int, int &gt; &amp;dst_idxs, const std::pair&lt; int, int &gt; &amp;src_idxs, int dst_ln)</arglist>
    </member>
    <member kind="function" protection="protected">
      <type>void</type>
      <name>xeqScheduleRestriction</name>
      <anchorfile>class_i_b_a_m_r_1_1_staggered_stokes_f_a_c_preconditioner_strategy.html</anchorfile>
      <anchor>a7a1fd1c33943ccb0eb4280d5e448648b</anchor>
      <arglist>(const std::pair&lt; int, int &gt; &amp;dst_idxs, const std::pair&lt; int, int &gt; &amp;src_idxs, int dst_ln)</arglist>
    </member>
    <member kind="function" protection="protected">
      <type>void</type>
      <name>xeqScheduleGhostFillNoCoarse</name>
      <anchorfile>class_i_b_a_m_r_1_1_staggered_stokes_f_a_c_preconditioner_strategy.html</anchorfile>
      <anchor>a9caef0cd2220d77127c1c27001bebd35</anchor>
      <arglist>(const std::pair&lt; int, int &gt; &amp;dst_idxs, int dst_ln)</arglist>
    </member>
    <member kind="function" protection="protected">
      <type>void</type>
      <name>xeqScheduleDataSynch</name>
      <anchorfile>class_i_b_a_m_r_1_1_staggered_stokes_f_a_c_preconditioner_strategy.html</anchorfile>
      <anchor>a07eeb5ecd09ed2714a4174e4b0b1e974</anchor>
      <arglist>(int dst_idx, int dst_ln)</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>IBAMR::StaggeredStokesFACPreconditioner</name>
    <filename>class_i_b_a_m_r_1_1_staggered_stokes_f_a_c_preconditioner.html</filename>
    <base>IBTK::FACPreconditioner</base>
    <base>IBAMR::StaggeredStokesSolver</base>
    <member kind="function">
      <type></type>
      <name>StaggeredStokesFACPreconditioner</name>
      <anchorfile>class_i_b_a_m_r_1_1_staggered_stokes_f_a_c_preconditioner.html</anchorfile>
      <anchor>ad5190391da63213b82c9390faa8fd2fb</anchor>
      <arglist>(const std::string &amp;object_name, SAMRAI::tbox::Pointer&lt; IBTK::FACPreconditionerStrategy &gt; fac_strategy, SAMRAI::tbox::Pointer&lt; SAMRAI::tbox::Database &gt; input_db, const std::string &amp;default_options_prefix)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~StaggeredStokesFACPreconditioner</name>
      <anchorfile>class_i_b_a_m_r_1_1_staggered_stokes_f_a_c_preconditioner.html</anchorfile>
      <anchor>a246e7859aed98afd5a64967ea0fd0b68</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setVelocityPoissonSpecifications</name>
      <anchorfile>class_i_b_a_m_r_1_1_staggered_stokes_f_a_c_preconditioner.html</anchorfile>
      <anchor>accaebaf7ce263f6efbbf9df64c8879e4</anchor>
      <arglist>(const SAMRAI::solv::PoissonSpecifications &amp;U_problem_coefs)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setPhysicalBcCoefs</name>
      <anchorfile>class_i_b_a_m_r_1_1_staggered_stokes_f_a_c_preconditioner.html</anchorfile>
      <anchor>abf28c52340f41b88f7c0d0572eb946ca</anchor>
      <arglist>(const std::vector&lt; SAMRAI::solv::RobinBcCoefStrategy&lt; NDIM &gt; * &gt; &amp;U_bc_coefs, SAMRAI::solv::RobinBcCoefStrategy&lt; NDIM &gt; *P_bc_coef)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setPhysicalBoundaryHelper</name>
      <anchorfile>class_i_b_a_m_r_1_1_staggered_stokes_f_a_c_preconditioner.html</anchorfile>
      <anchor>ac4537e1b0efa9b26a65c08ba9f5c63ac</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; StaggeredStokesPhysicalBoundaryHelper &gt; bc_helper)</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>IBAMR::StaggeredStokesFACPreconditionerStrategy</name>
    <filename>class_i_b_a_m_r_1_1_staggered_stokes_f_a_c_preconditioner_strategy.html</filename>
    <base>IBTK::FACPreconditionerStrategy</base>
    <member kind="function" protection="protected" virtualness="pure">
      <type>virtual void</type>
      <name>initializeOperatorStateSpecialized</name>
      <anchorfile>class_i_b_a_m_r_1_1_staggered_stokes_f_a_c_preconditioner_strategy.html</anchorfile>
      <anchor>a3dbe48e2ecf558945f8bc35603b670a6</anchor>
      <arglist>(const SAMRAI::solv::SAMRAIVectorReal&lt; NDIM, double &gt; &amp;solution, const SAMRAI::solv::SAMRAIVectorReal&lt; NDIM, double &gt; &amp;rhs, int coarsest_reset_ln, int finest_reset_ln)=0</arglist>
    </member>
    <member kind="function" protection="protected" virtualness="pure">
      <type>virtual void</type>
      <name>deallocateOperatorStateSpecialized</name>
      <anchorfile>class_i_b_a_m_r_1_1_staggered_stokes_f_a_c_preconditioner_strategy.html</anchorfile>
      <anchor>aba2c1f81af2223378240d9c5fac2ef1e</anchor>
      <arglist>(int coarsest_reset_ln, int finest_reset_ln)=0</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>IBAMR::StaggeredStokesOpenBoundaryStabilizer</name>
    <filename>class_i_b_a_m_r_1_1_staggered_stokes_open_boundary_stabilizer.html</filename>
    <base>IBTK::CartGridFunction</base>
    <member kind="function">
      <type></type>
      <name>StaggeredStokesOpenBoundaryStabilizer</name>
      <anchorfile>class_i_b_a_m_r_1_1_staggered_stokes_open_boundary_stabilizer.html</anchorfile>
      <anchor>a0189d02708a2fd85408a5f78ce883981</anchor>
      <arglist>(const std::string &amp;object_name, SAMRAI::tbox::Pointer&lt; SAMRAI::tbox::Database &gt; input_db, const INSHierarchyIntegrator *fluid_solver, SAMRAI::tbox::Pointer&lt; SAMRAI::geom::CartesianGridGeometry&lt; NDIM &gt; &gt; grid_geometry)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~StaggeredStokesOpenBoundaryStabilizer</name>
      <anchorfile>class_i_b_a_m_r_1_1_staggered_stokes_open_boundary_stabilizer.html</anchorfile>
      <anchor>af1cc1c7c419810358a8b7722e73d13d9</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>isTimeDependent</name>
      <anchorfile>class_i_b_a_m_r_1_1_staggered_stokes_open_boundary_stabilizer.html</anchorfile>
      <anchor>a31b9bc0cef9ad779cd76bdbf11b3a5a4</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setDataOnPatch</name>
      <anchorfile>class_i_b_a_m_r_1_1_staggered_stokes_open_boundary_stabilizer.html</anchorfile>
      <anchor>a8e9d7c50b0ca298306b071bbe5a8b713</anchor>
      <arglist>(int data_idx, SAMRAI::tbox::Pointer&lt; SAMRAI::hier::Variable&lt; NDIM &gt; &gt; var, SAMRAI::tbox::Pointer&lt; SAMRAI::hier::Patch&lt; NDIM &gt; &gt; patch, double data_time, bool initial_time=false, SAMRAI::tbox::Pointer&lt; SAMRAI::hier::PatchLevel&lt; NDIM &gt; &gt; level=SAMRAI::tbox::Pointer&lt; SAMRAI::hier::PatchLevel&lt; NDIM &gt; &gt;(NULL))</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>IBAMR::StaggeredStokesOperator</name>
    <filename>class_i_b_a_m_r_1_1_staggered_stokes_operator.html</filename>
    <base>IBTK::LinearOperator</base>
    <member kind="function">
      <type></type>
      <name>StaggeredStokesOperator</name>
      <anchorfile>class_i_b_a_m_r_1_1_staggered_stokes_operator.html</anchorfile>
      <anchor>ab824dd4e6c9944ff8e49f383ae0a49e8</anchor>
      <arglist>(const std::string &amp;object_name, bool homogeneous_bc=true)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~StaggeredStokesOperator</name>
      <anchorfile>class_i_b_a_m_r_1_1_staggered_stokes_operator.html</anchorfile>
      <anchor>a575b679025ff4b03030bdf58244195ee</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>setVelocityPoissonSpecifications</name>
      <anchorfile>class_i_b_a_m_r_1_1_staggered_stokes_operator.html</anchorfile>
      <anchor>ab3fe06329154ec987da2cdd9a2a5825e</anchor>
      <arglist>(const SAMRAI::solv::PoissonSpecifications &amp;U_problem_coefs)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual const SAMRAI::solv::PoissonSpecifications &amp;</type>
      <name>getVelocityPoissonSpecifications</name>
      <anchorfile>class_i_b_a_m_r_1_1_staggered_stokes_operator.html</anchorfile>
      <anchor>a01b2f45060a0e60091f232dc79f0e6f7</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>setPhysicalBcCoefs</name>
      <anchorfile>class_i_b_a_m_r_1_1_staggered_stokes_operator.html</anchorfile>
      <anchor>a7efed7a342bdada07125fc0efad106f7</anchor>
      <arglist>(const std::vector&lt; SAMRAI::solv::RobinBcCoefStrategy&lt; NDIM &gt; * &gt; &amp;U_bc_coefs, SAMRAI::solv::RobinBcCoefStrategy&lt; NDIM &gt; *P_bc_coef)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>setPhysicalBoundaryHelper</name>
      <anchorfile>class_i_b_a_m_r_1_1_staggered_stokes_operator.html</anchorfile>
      <anchor>a49bae721747d412481b66cca001be31a</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; StaggeredStokesPhysicalBoundaryHelper &gt; bc_helper)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>apply</name>
      <anchorfile>class_i_b_a_m_r_1_1_staggered_stokes_operator.html</anchorfile>
      <anchor>aac6e7d31a586cb59190f7053a5eecf3d</anchor>
      <arglist>(SAMRAI::solv::SAMRAIVectorReal&lt; NDIM, double &gt; &amp;x, SAMRAI::solv::SAMRAIVectorReal&lt; NDIM, double &gt; &amp;y)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>initializeOperatorState</name>
      <anchorfile>class_i_b_a_m_r_1_1_staggered_stokes_operator.html</anchorfile>
      <anchor>a3d90a8c2d0e841cdbd2f992cf5ae7e86</anchor>
      <arglist>(const SAMRAI::solv::SAMRAIVectorReal&lt; NDIM, double &gt; &amp;in, const SAMRAI::solv::SAMRAIVectorReal&lt; NDIM, double &gt; &amp;out)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>deallocateOperatorState</name>
      <anchorfile>class_i_b_a_m_r_1_1_staggered_stokes_operator.html</anchorfile>
      <anchor>aefc90a3e71f1d8144a36878b0b155591</anchor>
      <arglist>()</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>IBAMR::StaggeredStokesPETScLevelSolver</name>
    <filename>class_i_b_a_m_r_1_1_staggered_stokes_p_e_t_sc_level_solver.html</filename>
    <base>IBTK::PETScLevelSolver</base>
    <base>IBAMR::StaggeredStokesSolver</base>
    <member kind="function">
      <type></type>
      <name>StaggeredStokesPETScLevelSolver</name>
      <anchorfile>class_i_b_a_m_r_1_1_staggered_stokes_p_e_t_sc_level_solver.html</anchorfile>
      <anchor>aeb96ec4b900573512bdc66a77bb1c7d5</anchor>
      <arglist>(const std::string &amp;object_name, SAMRAI::tbox::Pointer&lt; SAMRAI::tbox::Database &gt; input_db, const std::string &amp;default_options_prefix)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~StaggeredStokesPETScLevelSolver</name>
      <anchorfile>class_i_b_a_m_r_1_1_staggered_stokes_p_e_t_sc_level_solver.html</anchorfile>
      <anchor>abed3268efe0155ba5d50d72c9241258e</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>setVelocityPoissonSpecifications</name>
      <anchorfile>class_i_b_a_m_r_1_1_staggered_stokes_solver.html</anchorfile>
      <anchor>a7a008276dfce88d0b19fefb1c9b90127</anchor>
      <arglist>(const SAMRAI::solv::PoissonSpecifications &amp;U_problem_coefs)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>setPhysicalBcCoefs</name>
      <anchorfile>class_i_b_a_m_r_1_1_staggered_stokes_solver.html</anchorfile>
      <anchor>aba34e52fc718270643d3f20791d1e107</anchor>
      <arglist>(const std::vector&lt; SAMRAI::solv::RobinBcCoefStrategy&lt; NDIM &gt; * &gt; &amp;U_bc_coefs, SAMRAI::solv::RobinBcCoefStrategy&lt; NDIM &gt; *P_bc_coef)</arglist>
    </member>
    <member kind="function" static="yes">
      <type>static SAMRAI::tbox::Pointer&lt; StaggeredStokesSolver &gt;</type>
      <name>allocate_solver</name>
      <anchorfile>class_i_b_a_m_r_1_1_staggered_stokes_p_e_t_sc_level_solver.html</anchorfile>
      <anchor>a97b57da466a3d6545c3d884ff95e182d</anchor>
      <arglist>(const std::string &amp;object_name, SAMRAI::tbox::Pointer&lt; SAMRAI::tbox::Database &gt; input_db, const std::string &amp;default_options_prefix)</arglist>
    </member>
    <member kind="function" protection="protected">
      <type>void</type>
      <name>initializeSolverStateSpecialized</name>
      <anchorfile>class_i_b_a_m_r_1_1_staggered_stokes_p_e_t_sc_level_solver.html</anchorfile>
      <anchor>a9589f5a27c3fea9c614bfa3140b6c261</anchor>
      <arglist>(const SAMRAI::solv::SAMRAIVectorReal&lt; NDIM, double &gt; &amp;x, const SAMRAI::solv::SAMRAIVectorReal&lt; NDIM, double &gt; &amp;b)</arglist>
    </member>
    <member kind="function" protection="protected">
      <type>void</type>
      <name>deallocateSolverStateSpecialized</name>
      <anchorfile>class_i_b_a_m_r_1_1_staggered_stokes_p_e_t_sc_level_solver.html</anchorfile>
      <anchor>ab9283d717b7f3110d0c5fef12d176de9</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" protection="protected">
      <type>void</type>
      <name>copyToPETScVec</name>
      <anchorfile>class_i_b_a_m_r_1_1_staggered_stokes_p_e_t_sc_level_solver.html</anchorfile>
      <anchor>a7f39206cb2bc1d1f2ec5e8fc252f130d</anchor>
      <arglist>(Vec &amp;petsc_x, SAMRAI::solv::SAMRAIVectorReal&lt; NDIM, double &gt; &amp;x, SAMRAI::tbox::Pointer&lt; SAMRAI::hier::PatchLevel&lt; NDIM &gt; &gt; patch_level)</arglist>
    </member>
    <member kind="function" protection="protected">
      <type>void</type>
      <name>copyFromPETScVec</name>
      <anchorfile>class_i_b_a_m_r_1_1_staggered_stokes_p_e_t_sc_level_solver.html</anchorfile>
      <anchor>a9a455abcb34d37d935e6a2a8635c67c5</anchor>
      <arglist>(Vec &amp;petsc_x, SAMRAI::solv::SAMRAIVectorReal&lt; NDIM, double &gt; &amp;x, SAMRAI::tbox::Pointer&lt; SAMRAI::hier::PatchLevel&lt; NDIM &gt; &gt; patch_level)</arglist>
    </member>
    <member kind="function" protection="protected">
      <type>void</type>
      <name>setupKSPVecs</name>
      <anchorfile>class_i_b_a_m_r_1_1_staggered_stokes_p_e_t_sc_level_solver.html</anchorfile>
      <anchor>aa402059d6e41bd86845fff0c6cc46b59</anchor>
      <arglist>(Vec &amp;petsc_x, Vec &amp;petsc_b, SAMRAI::solv::SAMRAIVectorReal&lt; NDIM, double &gt; &amp;x, SAMRAI::solv::SAMRAIVectorReal&lt; NDIM, double &gt; &amp;b, SAMRAI::tbox::Pointer&lt; SAMRAI::hier::PatchLevel&lt; NDIM &gt; &gt; patch_level)</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>IBAMR::StaggeredStokesPETScMatUtilities</name>
    <filename>class_i_b_a_m_r_1_1_staggered_stokes_p_e_t_sc_mat_utilities.html</filename>
    <member kind="function" static="yes">
      <type>static void</type>
      <name>constructPatchLevelMACStokesOp</name>
      <anchorfile>class_i_b_a_m_r_1_1_staggered_stokes_p_e_t_sc_mat_utilities.html</anchorfile>
      <anchor>adb227c8d1665d487687b9c2bc57d824e</anchor>
      <arglist>(Mat &amp;mat, const SAMRAI::solv::PoissonSpecifications &amp;u_problem_coefs, const std::vector&lt; SAMRAI::solv::RobinBcCoefStrategy&lt; NDIM &gt; * &gt; &amp;u_bc_coefs, double data_time, const std::vector&lt; int &gt; &amp;num_dofs_per_proc, int u_dof_index_idx, int p_dof_index_idx, SAMRAI::tbox::Pointer&lt; SAMRAI::hier::PatchLevel&lt; NDIM &gt; &gt; patch_level)</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>IBAMR::StaggeredStokesPETScVecUtilities</name>
    <filename>class_i_b_a_m_r_1_1_staggered_stokes_p_e_t_sc_vec_utilities.html</filename>
    <member kind="function" static="yes">
      <type>static void</type>
      <name>copyToPatchLevelVec</name>
      <anchorfile>class_i_b_a_m_r_1_1_staggered_stokes_p_e_t_sc_vec_utilities.html</anchorfile>
      <anchor>a090f4cad689cea4a393b22db4ed3061b</anchor>
      <arglist>(Vec &amp;vec, int u_data_idx, int u_dof_index_idx, int p_data_idx, int p_dof_index_idx, SAMRAI::tbox::Pointer&lt; SAMRAI::hier::PatchLevel&lt; NDIM &gt; &gt; patch_level)</arglist>
    </member>
    <member kind="function" static="yes">
      <type>static void</type>
      <name>copyFromPatchLevelVec</name>
      <anchorfile>class_i_b_a_m_r_1_1_staggered_stokes_p_e_t_sc_vec_utilities.html</anchorfile>
      <anchor>a426d02e96893026f883d7c7f00897bae</anchor>
      <arglist>(Vec &amp;vec, int u_data_idx, int u_dof_index_idx, int p_data_idx, int p_dof_index_idx, SAMRAI::tbox::Pointer&lt; SAMRAI::hier::PatchLevel&lt; NDIM &gt; &gt; patch_level, SAMRAI::tbox::Pointer&lt; SAMRAI::xfer::RefineSchedule&lt; NDIM &gt; &gt; data_synch_sched, SAMRAI::tbox::Pointer&lt; SAMRAI::xfer::RefineSchedule&lt; NDIM &gt; &gt; ghost_fill_sched)</arglist>
    </member>
    <member kind="function" static="yes">
      <type>static SAMRAI::tbox::Pointer&lt; SAMRAI::xfer::RefineSchedule&lt; NDIM &gt; &gt;</type>
      <name>constructDataSynchSchedule</name>
      <anchorfile>class_i_b_a_m_r_1_1_staggered_stokes_p_e_t_sc_vec_utilities.html</anchorfile>
      <anchor>af8fbee0cafeaa5473471a1786e150081</anchor>
      <arglist>(int u_data_idx, int p_data_idx, SAMRAI::tbox::Pointer&lt; SAMRAI::hier::PatchLevel&lt; NDIM &gt; &gt; patch_level)</arglist>
    </member>
    <member kind="function" static="yes">
      <type>static SAMRAI::tbox::Pointer&lt; SAMRAI::xfer::RefineSchedule&lt; NDIM &gt; &gt;</type>
      <name>constructGhostFillSchedule</name>
      <anchorfile>class_i_b_a_m_r_1_1_staggered_stokes_p_e_t_sc_vec_utilities.html</anchorfile>
      <anchor>ad64e16884a5227309d9ce0715123c7b9</anchor>
      <arglist>(int u_data_idx, int p_data_idx, SAMRAI::tbox::Pointer&lt; SAMRAI::hier::PatchLevel&lt; NDIM &gt; &gt; patch_level)</arglist>
    </member>
    <member kind="function" static="yes">
      <type>static void</type>
      <name>constructPatchLevelDOFIndices</name>
      <anchorfile>class_i_b_a_m_r_1_1_staggered_stokes_p_e_t_sc_vec_utilities.html</anchorfile>
      <anchor>a5d1c5bd55e5cd2d96af7f846d1fa949e</anchor>
      <arglist>(std::vector&lt; int &gt; &amp;num_dofs_per_proc, int u_dof_index_idx, int p_dof_index_idx, SAMRAI::tbox::Pointer&lt; SAMRAI::hier::PatchLevel&lt; NDIM &gt; &gt; patch_level)</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>IBAMR::StaggeredStokesPhysicalBoundaryHelper</name>
    <filename>class_i_b_a_m_r_1_1_staggered_stokes_physical_boundary_helper.html</filename>
    <base>IBTK::StaggeredPhysicalBoundaryHelper</base>
    <member kind="function">
      <type></type>
      <name>StaggeredStokesPhysicalBoundaryHelper</name>
      <anchorfile>class_i_b_a_m_r_1_1_staggered_stokes_physical_boundary_helper.html</anchorfile>
      <anchor>a1a34d0b107485b7d92ea941e4ec0da5e</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~StaggeredStokesPhysicalBoundaryHelper</name>
      <anchorfile>class_i_b_a_m_r_1_1_staggered_stokes_physical_boundary_helper.html</anchorfile>
      <anchor>a3d477662cf3b23a11bc7d35ad33cced6</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>enforceNormalVelocityBoundaryConditions</name>
      <anchorfile>class_i_b_a_m_r_1_1_staggered_stokes_physical_boundary_helper.html</anchorfile>
      <anchor>aa9ab1b40627d135a28d96598344fd96b</anchor>
      <arglist>(int u_data_idx, int p_data_idx, const std::vector&lt; SAMRAI::solv::RobinBcCoefStrategy&lt; NDIM &gt; * &gt; &amp;u_bc_coefs, double fill_time, bool homogeneous_bc, int coarsest_ln=-1, int finest_ln=-1) const </arglist>
    </member>
    <member kind="function" static="yes">
      <type>static void</type>
      <name>setupBcCoefObjects</name>
      <anchorfile>class_i_b_a_m_r_1_1_staggered_stokes_physical_boundary_helper.html</anchorfile>
      <anchor>a55c89f3c51d3cf986af83ae5a66cbe1d</anchor>
      <arglist>(const std::vector&lt; SAMRAI::solv::RobinBcCoefStrategy&lt; NDIM &gt; * &gt; &amp;u_bc_coefs, SAMRAI::solv::RobinBcCoefStrategy&lt; NDIM &gt; *p_bc_coef, int u_target_data_idx, int p_target_data_idx, bool homogeneous_bc)</arglist>
    </member>
    <member kind="function" static="yes">
      <type>static void</type>
      <name>resetBcCoefObjects</name>
      <anchorfile>class_i_b_a_m_r_1_1_staggered_stokes_physical_boundary_helper.html</anchorfile>
      <anchor>a11a4d38f08c4444c0ff2578e22be8155</anchor>
      <arglist>(const std::vector&lt; SAMRAI::solv::RobinBcCoefStrategy&lt; NDIM &gt; * &gt; &amp;u_bc_coefs, SAMRAI::solv::RobinBcCoefStrategy&lt; NDIM &gt; *p_bc_coef)</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>IBAMR::StaggeredStokesProjectionPreconditioner</name>
    <filename>class_i_b_a_m_r_1_1_staggered_stokes_projection_preconditioner.html</filename>
    <base>IBAMR::StaggeredStokesBlockPreconditioner</base>
    <member kind="function">
      <type></type>
      <name>StaggeredStokesProjectionPreconditioner</name>
      <anchorfile>class_i_b_a_m_r_1_1_staggered_stokes_projection_preconditioner.html</anchorfile>
      <anchor>a7e93144d2d19e38bc077aec2ce57357f</anchor>
      <arglist>(const std::string &amp;object_name, SAMRAI::tbox::Pointer&lt; SAMRAI::tbox::Database &gt; input_db, const std::string &amp;default_options_prefix)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~StaggeredStokesProjectionPreconditioner</name>
      <anchorfile>class_i_b_a_m_r_1_1_staggered_stokes_projection_preconditioner.html</anchorfile>
      <anchor>a1dd7115b232a3f1e79b57d61b4476ca2</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>solveSystem</name>
      <anchorfile>class_i_b_a_m_r_1_1_staggered_stokes_projection_preconditioner.html</anchorfile>
      <anchor>ac37eda37ce49160a21f9fc08ebf81c62</anchor>
      <arglist>(SAMRAI::solv::SAMRAIVectorReal&lt; NDIM, double &gt; &amp;x, SAMRAI::solv::SAMRAIVectorReal&lt; NDIM, double &gt; &amp;b)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>initializeSolverState</name>
      <anchorfile>class_i_b_a_m_r_1_1_staggered_stokes_projection_preconditioner.html</anchorfile>
      <anchor>a114bbd528c76c82777c5099113d38fce</anchor>
      <arglist>(const SAMRAI::solv::SAMRAIVectorReal&lt; NDIM, double &gt; &amp;x, const SAMRAI::solv::SAMRAIVectorReal&lt; NDIM, double &gt; &amp;b)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>deallocateSolverState</name>
      <anchorfile>class_i_b_a_m_r_1_1_staggered_stokes_projection_preconditioner.html</anchorfile>
      <anchor>a6d3ca6a4b86fe33e3ac58167670b01df</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setInitialGuessNonzero</name>
      <anchorfile>class_i_b_a_m_r_1_1_staggered_stokes_projection_preconditioner.html</anchorfile>
      <anchor>a9bb9f64acacf162d5564702f2dd64805</anchor>
      <arglist>(bool initial_guess_nonzero=true)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setMaxIterations</name>
      <anchorfile>class_i_b_a_m_r_1_1_staggered_stokes_projection_preconditioner.html</anchorfile>
      <anchor>ae43c3367caa54fcb9bad3f778e170799</anchor>
      <arglist>(int max_iterations)</arglist>
    </member>
    <member kind="function" static="yes">
      <type>static SAMRAI::tbox::Pointer&lt; StaggeredStokesSolver &gt;</type>
      <name>allocate_solver</name>
      <anchorfile>class_i_b_a_m_r_1_1_staggered_stokes_projection_preconditioner.html</anchorfile>
      <anchor>a4e90498d209c7adae5ed79e4d5bae415</anchor>
      <arglist>(const std::string &amp;object_name, SAMRAI::tbox::Pointer&lt; SAMRAI::tbox::Database &gt; input_db, const std::string &amp;default_options_prefix)</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>IBAMR::StaggeredStokesSolver</name>
    <filename>class_i_b_a_m_r_1_1_staggered_stokes_solver.html</filename>
    <base virtualness="virtual">IBTK::GeneralSolver</base>
  </compound>
  <compound kind="class">
    <name>IBAMR::StaggeredStokesSolverManager</name>
    <filename>class_i_b_a_m_r_1_1_staggered_stokes_solver_manager.html</filename>
    <member kind="typedef">
      <type>SAMRAI::tbox::Pointer&lt; StaggeredStokesSolver &gt;(*</type>
      <name>SolverMaker</name>
      <anchorfile>class_i_b_a_m_r_1_1_staggered_stokes_solver_manager.html</anchorfile>
      <anchor>ad1e580ec0f659a69aa1e6a58594ea09e</anchor>
      <arglist>)(const std::string &amp;solver_object_name, SAMRAI::tbox::Pointer&lt; SAMRAI::tbox::Database &gt; solver_input_db, const std::string &amp;solver_default_options_prefix)</arglist>
    </member>
    <member kind="function">
      <type>SAMRAI::tbox::Pointer&lt; StaggeredStokesSolver &gt;</type>
      <name>allocateSolver</name>
      <anchorfile>class_i_b_a_m_r_1_1_staggered_stokes_solver_manager.html</anchorfile>
      <anchor>a76378a194adf9b781021607b0b249dd6</anchor>
      <arglist>(const std::string &amp;solver_type, const std::string &amp;solver_object_name, SAMRAI::tbox::Pointer&lt; SAMRAI::tbox::Database &gt; solver_input_db, const std::string &amp;solver_default_options_prefix) const </arglist>
    </member>
    <member kind="function">
      <type>SAMRAI::tbox::Pointer&lt; StaggeredStokesSolver &gt;</type>
      <name>allocateSolver</name>
      <anchorfile>class_i_b_a_m_r_1_1_staggered_stokes_solver_manager.html</anchorfile>
      <anchor>a4d1dfbc139bb43513a4e713b7b2763ec</anchor>
      <arglist>(const std::string &amp;solver_type, const std::string &amp;solver_object_name, SAMRAI::tbox::Pointer&lt; SAMRAI::tbox::Database &gt; solver_input_db, const std::string &amp;solver_default_options_prefix, const std::string &amp;precond_type, const std::string &amp;precond_object_name, SAMRAI::tbox::Pointer&lt; SAMRAI::tbox::Database &gt; precond_input_db, const std::string &amp;precond_default_options_prefix) const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>registerSolverFactoryFunction</name>
      <anchorfile>class_i_b_a_m_r_1_1_staggered_stokes_solver_manager.html</anchorfile>
      <anchor>a67714278f6a6d72c09f9548ef21e52d4</anchor>
      <arglist>(const std::string &amp;solver_type, SolverMaker solver_maker)</arglist>
    </member>
    <member kind="function" static="yes">
      <type>static StaggeredStokesSolverManager *</type>
      <name>getManager</name>
      <anchorfile>class_i_b_a_m_r_1_1_staggered_stokes_solver_manager.html</anchorfile>
      <anchor>af3e7e739d41a5d65fd75141a33141f00</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" static="yes">
      <type>static void</type>
      <name>freeManager</name>
      <anchorfile>class_i_b_a_m_r_1_1_staggered_stokes_solver_manager.html</anchorfile>
      <anchor>af26de0809ce78814e8b905c103def4ce</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="variable" static="yes">
      <type>static const std::string</type>
      <name>UNDEFINED</name>
      <anchorfile>class_i_b_a_m_r_1_1_staggered_stokes_solver_manager.html</anchorfile>
      <anchor>abbf622a4cc45f5397d643da13a5a800a</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" static="yes">
      <type>static const std::string</type>
      <name>DEFAULT_KRYLOV_SOLVER</name>
      <anchorfile>class_i_b_a_m_r_1_1_staggered_stokes_solver_manager.html</anchorfile>
      <anchor>aada187888c0688a606aa5ec5d7df7c44</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" static="yes">
      <type>static const std::string</type>
      <name>DEFAULT_BLOCK_PRECONDITIONER</name>
      <anchorfile>class_i_b_a_m_r_1_1_staggered_stokes_solver_manager.html</anchorfile>
      <anchor>a44ed14163261bd0885cdb634b8650650</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" static="yes">
      <type>static const std::string</type>
      <name>DEFAULT_FAC_PRECONDITIONER</name>
      <anchorfile>class_i_b_a_m_r_1_1_staggered_stokes_solver_manager.html</anchorfile>
      <anchor>a91d7920421fe0d61a4ca92803ed2b12c</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" static="yes">
      <type>static const std::string</type>
      <name>DEFAULT_LEVEL_SOLVER</name>
      <anchorfile>class_i_b_a_m_r_1_1_staggered_stokes_solver_manager.html</anchorfile>
      <anchor>a1192aef7ef1cb5951e3363277f55c225</anchor>
      <arglist></arglist>
    </member>
    <member kind="function" protection="protected">
      <type></type>
      <name>StaggeredStokesSolverManager</name>
      <anchorfile>class_i_b_a_m_r_1_1_staggered_stokes_solver_manager.html</anchorfile>
      <anchor>a39dd79e673f3c556ebcef07cc8cc588f</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" protection="protected">
      <type></type>
      <name>~StaggeredStokesSolverManager</name>
      <anchorfile>class_i_b_a_m_r_1_1_staggered_stokes_solver_manager.html</anchorfile>
      <anchor>abbbe7149e099a4b923d6ec1f9909d016</anchor>
      <arglist>()</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>IBAMR::StokesBcCoefStrategy</name>
    <filename>class_i_b_a_m_r_1_1_stokes_bc_coef_strategy.html</filename>
    <base>IBTK::ExtendedRobinBcCoefStrategy</base>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>setStokesSpecifications</name>
      <anchorfile>class_i_b_a_m_r_1_1_stokes_bc_coef_strategy.html</anchorfile>
      <anchor>a068575db447f8fdbd9bf28e3b78e2e9b</anchor>
      <arglist>(const StokesSpecifications *problem_coefs)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>setTargetVelocityPatchDataIndex</name>
      <anchorfile>class_i_b_a_m_r_1_1_stokes_bc_coef_strategy.html</anchorfile>
      <anchor>a122114aa6b6dabd663c5e0ba568e2fc6</anchor>
      <arglist>(int u_target_data_idx)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>clearTargetVelocityPatchDataIndex</name>
      <anchorfile>class_i_b_a_m_r_1_1_stokes_bc_coef_strategy.html</anchorfile>
      <anchor>a757844997bd0b40c7ef92790d4b285b8</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>setTargetPressurePatchDataIndex</name>
      <anchorfile>class_i_b_a_m_r_1_1_stokes_bc_coef_strategy.html</anchorfile>
      <anchor>a6b64e5a7954af10ea4f0336e9e2e3df9</anchor>
      <arglist>(int u_target_data_idx)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>clearTargetPressurePatchDataIndex</name>
      <anchorfile>class_i_b_a_m_r_1_1_stokes_bc_coef_strategy.html</anchorfile>
      <anchor>a999f270844f98c840f7018e50e68413c</anchor>
      <arglist>()</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>IBAMR::StokesSpecifications</name>
    <filename>class_i_b_a_m_r_1_1_stokes_specifications.html</filename>
    <member kind="function">
      <type></type>
      <name>StokesSpecifications</name>
      <anchorfile>class_i_b_a_m_r_1_1_stokes_specifications.html</anchorfile>
      <anchor>ac68cbdf5b4e2d5dc5985eedaccaea78a</anchor>
      <arglist>(double rho=0.0, double mu=0.0, double lambda=0.0)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>StokesSpecifications</name>
      <anchorfile>class_i_b_a_m_r_1_1_stokes_specifications.html</anchorfile>
      <anchor>adada8233bbfd4448ac8e193c8729bb09</anchor>
      <arglist>(const StokesSpecifications &amp;from)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~StokesSpecifications</name>
      <anchorfile>class_i_b_a_m_r_1_1_stokes_specifications.html</anchorfile>
      <anchor>a03b97935c6752f89024b37ce30c42d87</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>StokesSpecifications &amp;</type>
      <name>operator=</name>
      <anchorfile>class_i_b_a_m_r_1_1_stokes_specifications.html</anchorfile>
      <anchor>a4f48c6fe6b6288727febca6f70ecd48d</anchor>
      <arglist>(const StokesSpecifications &amp;that)</arglist>
    </member>
    <member kind="function">
      <type>double</type>
      <name>getRho</name>
      <anchorfile>class_i_b_a_m_r_1_1_stokes_specifications.html</anchorfile>
      <anchor>a67ee234284ad5aa05cc9d2e1f249fd5c</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setRho</name>
      <anchorfile>class_i_b_a_m_r_1_1_stokes_specifications.html</anchorfile>
      <anchor>a9dfc16011dbc02d7cb81caa78f9f92bd</anchor>
      <arglist>(double rho)</arglist>
    </member>
    <member kind="function">
      <type>double</type>
      <name>getMu</name>
      <anchorfile>class_i_b_a_m_r_1_1_stokes_specifications.html</anchorfile>
      <anchor>a85cf20bb42af5a4611e5f322594c3695</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setMu</name>
      <anchorfile>class_i_b_a_m_r_1_1_stokes_specifications.html</anchorfile>
      <anchor>ad1a891fa982b4d564fee4eab333baa47</anchor>
      <arglist>(double mu)</arglist>
    </member>
    <member kind="function">
      <type>double</type>
      <name>getLambda</name>
      <anchorfile>class_i_b_a_m_r_1_1_stokes_specifications.html</anchorfile>
      <anchor>a0421da55634de38db7d80e5227bb3fa1</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setLambda</name>
      <anchorfile>class_i_b_a_m_r_1_1_stokes_specifications.html</anchorfile>
      <anchor>a3f2c75eb4b035b97850cead2d0be9493</anchor>
      <arglist>(double lambda)</arglist>
    </member>
  </compound>
  <compound kind="namespace">
    <name>IBTK</name>
    <filename>namespace_i_b_t_k.html</filename>
  </compound>
  <compound kind="struct">
    <name>std::less&lt; SAMRAI::tbox::Pointer&lt; T &gt; &gt;</name>
    <filename>structstd_1_1less_3_01_s_a_m_r_a_i_1_1tbox_1_1_pointer_3_01_t_01_4_01_4.html</filename>
    <templarg></templarg>
  </compound>
  <compound kind="dir">
    <name>/Users/griffith/code/IBAMR/src/adv_diff</name>
    <path>/Users/griffith/code/IBAMR/src/adv_diff/</path>
    <filename>dir_1f0b91ccdeba06ce428a19b8ef65e6d0.html</filename>
    <file>AdvDiffCenteredConvectiveOperator.cpp</file>
    <file>AdvDiffConvectiveOperatorManager.cpp</file>
    <file>AdvDiffHierarchyIntegrator.cpp</file>
    <file>AdvDiffPhysicalBoundaryUtilities.cpp</file>
    <file>AdvDiffPPMConvectiveOperator.cpp</file>
    <file>AdvDiffPredictorCorrectorHierarchyIntegrator.cpp</file>
    <file>AdvDiffPredictorCorrectorHyperbolicPatchOps.cpp</file>
    <file>AdvDiffSemiImplicitHierarchyIntegrator.cpp</file>
    <file>AdvDiffStochasticForcing.cpp</file>
  </compound>
  <compound kind="dir">
    <name>/Users/griffith/code/IBAMR/src/advect</name>
    <path>/Users/griffith/code/IBAMR/src/advect/</path>
    <filename>dir_0023837800241f8a50a9944050ae0c4a.html</filename>
    <file>AdvectorExplicitPredictorPatchOps.cpp</file>
    <file>AdvectorPredictorCorrectorHyperbolicPatchOps.cpp</file>
  </compound>
  <compound kind="dir">
    <name>/Users/griffith/code/IBAMR/src/IB</name>
    <path>/Users/griffith/code/IBAMR/src/IB/</path>
    <filename>dir_a7f6aa32002053ae7b20aefbbc6b429a.html</filename>
    <file>GeneralizedIBMethod.cpp</file>
    <file>IBAnchorPointSpec.cpp</file>
    <file>IBAnchorPointSpecFactory.cpp</file>
    <file>IBBeamForceSpec.cpp</file>
    <file>IBBeamForceSpecFactory.cpp</file>
    <file>IBEulerianForceFunction.cpp</file>
    <file>IBEulerianSourceFunction.cpp</file>
    <file>IBExplicitHierarchyIntegrator.cpp</file>
    <file>IBFECentroidPostProcessor.cpp</file>
    <file>IBFEMethod.cpp</file>
    <file>IBFEPatchRecoveryPostProcessor.cpp</file>
    <file>IBFEPostProcessor.cpp</file>
    <file>IBHierarchyIntegrator.cpp</file>
    <file>IBImplicitStaggeredHierarchyIntegrator.cpp</file>
    <file>IBImplicitStrategy.cpp</file>
    <file>IBInstrumentationSpec.cpp</file>
    <file>IBInstrumentationSpecFactory.cpp</file>
    <file>IBInstrumentPanel.cpp</file>
    <file>IBKirchhoffRodForceGen.cpp</file>
    <file>IBLagrangianForceStrategy.cpp</file>
    <file>IBLagrangianForceStrategySet.cpp</file>
    <file>IBLagrangianSourceStrategy.cpp</file>
    <file>IBMethod.cpp</file>
    <file>IBMethodPostProcessStrategy.cpp</file>
    <file>IBRodForceSpec.cpp</file>
    <file>IBRodForceSpecFactory.cpp</file>
    <file>IBSourceSpec.cpp</file>
    <file>IBSourceSpecFactory.cpp</file>
    <file>IBSpringForceSpec.cpp</file>
    <file>IBSpringForceSpecFactory.cpp</file>
    <file>IBStandardForceGen.cpp</file>
    <file>IBStandardInitializer.cpp</file>
    <file>IBStandardSourceGen.cpp</file>
    <file>IBStrategy.cpp</file>
    <file>IBStrategySet.cpp</file>
    <file>IBTargetPointForceSpec.cpp</file>
    <file>IBTargetPointForceSpecFactory.cpp</file>
    <file>IMPInitializer.cpp</file>
    <file>IMPMethod.cpp</file>
    <file>MaterialPointSpec.cpp</file>
    <file>MaterialPointSpecFactory.cpp</file>
    <file>PenaltyIBMethod.cpp</file>
  </compound>
  <compound kind="dir">
    <name>/Users/griffith/code/IBAMR/include/ibamr</name>
    <path>/Users/griffith/code/IBAMR/include/ibamr/</path>
    <filename>dir_9e978f3cc14d40f2c17d3563240b4162.html</filename>
    <file>AdvDiffCenteredConvectiveOperator.h</file>
    <file>AdvDiffConvectiveOperatorManager.h</file>
    <file>AdvDiffHierarchyIntegrator.h</file>
    <file>AdvDiffPhysicalBoundaryUtilities.h</file>
    <file>AdvDiffPPMConvectiveOperator.h</file>
    <file>AdvDiffPredictorCorrectorHierarchyIntegrator.h</file>
    <file>AdvDiffPredictorCorrectorHyperbolicPatchOps.h</file>
    <file>AdvDiffSemiImplicitHierarchyIntegrator.h</file>
    <file>AdvDiffStochasticForcing.h</file>
    <file>AdvectorExplicitPredictorPatchOps.h</file>
    <file>AdvectorPredictorCorrectorHyperbolicPatchOps.h</file>
    <file>app_namespaces.h</file>
    <file>ConvectiveOperator.h</file>
    <file>GeneralizedIBMethod.h</file>
    <file>ibamr.h</file>
    <file>ibamr_enums.h</file>
    <file>ibamr_utilities.h</file>
    <file>IBAnchorPointSpec-inl.h</file>
    <file>IBAnchorPointSpec.h</file>
    <file>IBBeamForceSpec-inl.h</file>
    <file>IBBeamForceSpec.h</file>
    <file>IBExplicitHierarchyIntegrator.h</file>
    <file>IBFECentroidPostProcessor.h</file>
    <file>IBFEMethod.h</file>
    <file>IBFEPatchRecoveryPostProcessor.h</file>
    <file>IBFEPostProcessor.h</file>
    <file>IBHierarchyIntegrator.h</file>
    <file>IBImplicitStaggeredHierarchyIntegrator.h</file>
    <file>IBImplicitStrategy.h</file>
    <file>IBInstrumentationSpec-inl.h</file>
    <file>IBInstrumentationSpec.h</file>
    <file>IBInstrumentPanel.h</file>
    <file>IBKirchhoffRodForceGen.h</file>
    <file>IBLagrangianForceStrategy.h</file>
    <file>IBLagrangianForceStrategySet.h</file>
    <file>IBLagrangianSourceStrategy.h</file>
    <file>IBMethod.h</file>
    <file>IBMethodPostProcessStrategy.h</file>
    <file>IBRodForceSpec-inl.h</file>
    <file>IBRodForceSpec.h</file>
    <file>IBSourceSpec-inl.h</file>
    <file>IBSourceSpec.h</file>
    <file>IBSpringForceFunctions.h</file>
    <file>IBSpringForceSpec-inl.h</file>
    <file>IBSpringForceSpec.h</file>
    <file>IBStandardForceGen.h</file>
    <file>IBStandardInitializer.h</file>
    <file>IBStandardSourceGen.h</file>
    <file>IBStrategy.h</file>
    <file>IBStrategySet.h</file>
    <file>IBTargetPointForceSpec-inl.h</file>
    <file>IBTargetPointForceSpec.h</file>
    <file>IMPInitializer.h</file>
    <file>IMPMethod.h</file>
    <file>INSCollocatedCenteredConvectiveOperator.h</file>
    <file>INSCollocatedConvectiveOperatorManager.h</file>
    <file>INSCollocatedHierarchyIntegrator.h</file>
    <file>INSCollocatedPPMConvectiveOperator.h</file>
    <file>INSCollocatedVelocityBcCoef.h</file>
    <file>INSHierarchyIntegrator.h</file>
    <file>INSIntermediateVelocityBcCoef.h</file>
    <file>INSProjectionBcCoef.h</file>
    <file>INSStaggeredCenteredConvectiveOperator.h</file>
    <file>INSStaggeredConvectiveOperatorManager.h</file>
    <file>INSStaggeredHierarchyIntegrator.h</file>
    <file>INSStaggeredPPMConvectiveOperator.h</file>
    <file>INSStaggeredPressureBcCoef.h</file>
    <file>INSStaggeredStabilizedPPMConvectiveOperator.h</file>
    <file>INSStaggeredStochasticForcing.h</file>
    <file>INSStaggeredUpwindConvectiveOperator.h</file>
    <file>INSStaggeredVelocityBcCoef.h</file>
    <file>KrylovLinearSolverStaggeredStokesSolverInterface.h</file>
    <file>MaterialPointSpec-inl.h</file>
    <file>MaterialPointSpec.h</file>
    <file>namespaces.h</file>
    <file>PenaltyIBMethod.h</file>
    <file>PETScKrylovStaggeredStokesSolver.h</file>
    <file>RNG.h</file>
    <file>SpongeLayerForceFunction.h</file>
    <file>StaggeredStokesBlockFactorizationPreconditioner.h</file>
    <file>StaggeredStokesBlockPreconditioner.h</file>
    <file>StaggeredStokesBoxRelaxationFACOperator.h</file>
    <file>StaggeredStokesFACPreconditioner.h</file>
    <file>StaggeredStokesFACPreconditionerStrategy.h</file>
    <file>StaggeredStokesOpenBoundaryStabilizer.h</file>
    <file>StaggeredStokesOperator.h</file>
    <file>StaggeredStokesPETScLevelSolver.h</file>
    <file>StaggeredStokesPETScMatUtilities.h</file>
    <file>StaggeredStokesPETScVecUtilities.h</file>
    <file>StaggeredStokesPhysicalBoundaryHelper.h</file>
    <file>StaggeredStokesProjectionPreconditioner.h</file>
    <file>StaggeredStokesSolver.h</file>
    <file>StaggeredStokesSolverManager.h</file>
    <file>StokesBcCoefStrategy.h</file>
    <file>StokesSpecifications.h</file>
  </compound>
  <compound kind="dir">
    <name>/Users/griffith/code/IBAMR/include</name>
    <path>/Users/griffith/code/IBAMR/include/</path>
    <filename>dir_d44c64559bbebec7f509842c48db8b23.html</filename>
    <dir>/Users/griffith/code/IBAMR/include/ibamr</dir>
  </compound>
  <compound kind="dir">
    <name>/Users/griffith/code/IBAMR/src/navier_stokes</name>
    <path>/Users/griffith/code/IBAMR/src/navier_stokes/</path>
    <filename>dir_0a45b23f7299b2a0eac0f5b1b08d0c04.html</filename>
    <file>INSCollocatedCenteredConvectiveOperator.cpp</file>
    <file>INSCollocatedConvectiveOperatorManager.cpp</file>
    <file>INSCollocatedHierarchyIntegrator.cpp</file>
    <file>INSCollocatedPPMConvectiveOperator.cpp</file>
    <file>INSCollocatedVelocityBcCoef.cpp</file>
    <file>INSHierarchyIntegrator.cpp</file>
    <file>INSIntermediateVelocityBcCoef.cpp</file>
    <file>INSProjectionBcCoef.cpp</file>
    <file>INSStaggeredCenteredConvectiveOperator.cpp</file>
    <file>INSStaggeredConvectiveOperatorManager.cpp</file>
    <file>INSStaggeredHierarchyIntegrator.cpp</file>
    <file>INSStaggeredPPMConvectiveOperator.cpp</file>
    <file>INSStaggeredPressureBcCoef.cpp</file>
    <file>INSStaggeredStabilizedPPMConvectiveOperator.cpp</file>
    <file>INSStaggeredStochasticForcing.cpp</file>
    <file>INSStaggeredUpwindConvectiveOperator.cpp</file>
    <file>INSStaggeredVelocityBcCoef.cpp</file>
    <file>KrylovLinearSolverStaggeredStokesSolverInterface.cpp</file>
    <file>PETScKrylovStaggeredStokesSolver.cpp</file>
    <file>SpongeLayerForceFunction.cpp</file>
    <file>StaggeredStokesBlockFactorizationPreconditioner.cpp</file>
    <file>StaggeredStokesBlockPreconditioner.cpp</file>
    <file>StaggeredStokesBoxRelaxationFACOperator.cpp</file>
    <file>StaggeredStokesFACPreconditioner.cpp</file>
    <file>StaggeredStokesFACPreconditionerStrategy.cpp</file>
    <file>StaggeredStokesOpenBoundaryStabilizer.cpp</file>
    <file>StaggeredStokesOperator.cpp</file>
    <file>StaggeredStokesPETScLevelSolver.cpp</file>
    <file>StaggeredStokesPETScMatUtilities.cpp</file>
    <file>StaggeredStokesPETScVecUtilities.cpp</file>
    <file>StaggeredStokesPhysicalBoundaryHelper.cpp</file>
    <file>StaggeredStokesProjectionPreconditioner.cpp</file>
    <file>StaggeredStokesSolver.cpp</file>
    <file>StaggeredStokesSolverManager.cpp</file>
    <file>StokesBcCoefStrategy.cpp</file>
    <file>StokesSpecifications.cpp</file>
  </compound>
  <compound kind="dir">
    <name>/Users/griffith/code/IBAMR/src</name>
    <path>/Users/griffith/code/IBAMR/src/</path>
    <filename>dir_68267d1309a1af8e8297ef4c3efbcdba.html</filename>
    <dir>/Users/griffith/code/IBAMR/src/adv_diff</dir>
    <dir>/Users/griffith/code/IBAMR/src/advect</dir>
    <dir>/Users/griffith/code/IBAMR/src/IB</dir>
    <dir>/Users/griffith/code/IBAMR/src/navier_stokes</dir>
    <dir>/Users/griffith/code/IBAMR/src/utilities</dir>
  </compound>
  <compound kind="dir">
    <name>/Users/griffith/code/IBAMR/src/utilities</name>
    <path>/Users/griffith/code/IBAMR/src/utilities/</path>
    <filename>dir_7b5d38f1875f1b693f62ca6a108a1129.html</filename>
    <file>ConvectiveOperator.cpp</file>
    <file>RNG.cpp</file>
  </compound>
  <compound kind="page">
    <name>index</name>
    <title></title>
    <filename>index</filename>
  </compound>
</tagfile>
