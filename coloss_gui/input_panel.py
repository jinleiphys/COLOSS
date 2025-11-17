"""
Input panel for COLOSS parameters with tabbed interface
"""

from PySide6.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QGroupBox, QLabel,
    QLineEdit, QSpinBox, QDoubleSpinBox, QCheckBox,
    QFormLayout, QPushButton, QTabWidget, QGridLayout, QFileDialog, QMessageBox
)
from PySide6.QtCore import Qt, Signal
from PySide6.QtGui import QFont, QCursor


class ClickableLabel(QLabel):
    """A clickable label that shows a message box on click"""
    def __init__(self, text, message, parent=None):
        super().__init__(text, parent)
        self.message = message
        self.setCursor(QCursor(Qt.PointingHandCursor))

    def mousePressEvent(self, event):
        if event.button() == Qt.LeftButton:
            QMessageBox.information(self, "Parameter Information", self.message)
        super().mousePressEvent(event)


class InputPanel(QWidget):
    """Panel for entering COLOSS input parameters"""

    input_changed = Signal()

    def __init__(self):
        super().__init__()
        self.pot_file_path = None  # Store path to external potential file
        self.init_ui()

    def create_label_with_info(self, text, tooltip):
        """Create a label with an info icon that shows tooltip on hover and click"""
        container = QWidget()
        layout = QHBoxLayout(container)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(5)

        # Main label
        label = QLabel(text)
        layout.addWidget(label)

        # Info icon with tooltip and click support
        info_icon = ClickableLabel("ⓘ", tooltip)
        info_icon.setStyleSheet("color: #3498db; font-size: 14px;")
        info_icon.setToolTip(tooltip)
        layout.addWidget(info_icon)

        layout.addStretch()
        return container

    def init_ui(self):
        """Initialize the input panel with tabbed interface"""
        main_layout = QVBoxLayout(self)
        main_layout.setContentsMargins(10, 10, 10, 10)

        # Create tab widget
        self.tabs = QTabWidget()

        # Create individual tabs
        self.tabs.addTab(self.create_general_tab(), "General")
        self.tabs.addTab(self.create_system_tab(), "System")
        self.tabs.addTab(self.create_potential_tab(), "Potential")
        self.tabs.addTab(self.create_nonlocal_tab(), "Nonlocal")

        main_layout.addWidget(self.tabs)

        # Buttons at bottom
        button_layout = QHBoxLayout()

        self.load_btn = QPushButton("Load from File")
        button_layout.addWidget(self.load_btn)

        button_layout.addStretch()

        self.generate_btn = QPushButton("Generate Input File")
        self.generate_btn.clicked.connect(self.generate_input)
        button_layout.addWidget(self.generate_btn)

        main_layout.addLayout(button_layout)

    def create_general_tab(self):
        """Create general parameters tab"""
        widget = QWidget()
        layout = QVBoxLayout(widget)

        # Title
        title = QLabel("General Parameters")
        title.setStyleSheet("font-size: 16px; font-weight: bold;")
        title.setAlignment(Qt.AlignCenter)
        layout.addWidget(title)

        # Two-column grid
        grid = QGridLayout()
        grid.setColumnStretch(0, 1)
        grid.setColumnStretch(1, 1)

        # Left column
        left_form = QFormLayout()

        self.ctheta = QSpinBox()
        self.ctheta.setRange(1, 10)
        self.ctheta.setValue(7)
        ctheta_label = self.create_label_with_info(
            "Rotation Angle (deg):",
            "Contour rotation angle θ (in degrees)\n"
            "Rotates integration contour: r → r·exp(iθ)\n"
            "Typical values: 3-7 degrees"
        )
        left_form.addRow(ctheta_label, self.ctheta)

        self.alpha = QDoubleSpinBox()
        self.alpha.setRange(-10, 10)
        self.alpha.setValue(0)
        alpha_label = self.create_label_with_info(
            "α parameter:",
            "α parameter in Generalized Laguerre polynomials L_n^(α)(x)\n\n"
            "Defines the weight function: w(x) = x^α · exp(-x)\n"
            "Used for Gauss-Laguerre quadrature on [0, ∞)\n\n"
            "Effects:\n"
            "  α = 0: Standard Laguerre (uniform weight decay)\n"
            "  α > 0: Emphasizes larger r (better for scattering)\n"
            "  α < 0: Emphasizes smaller r (better for bound states)\n\n"
            "Typical: α = 0 (standard Laguerre polynomials)"
        )
        left_form.addRow(alpha_label, self.alpha)

        self.nr = QSpinBox()
        self.nr.setRange(10, 1000)
        self.nr.setValue(60)
        self.nr.setToolTip("Number of radial basis functions")
        left_form.addRow("No. of basis functions:", self.nr)

        self.rmax = QDoubleSpinBox()
        self.rmax.setRange(1, 200)
        self.rmax.setValue(40)
        self.rmax.setToolTip("Maximum radial distance (fm)")
        self.rmax.valueChanged.connect(self.on_rmax_changed)
        left_form.addRow("Rmax (fm):", self.rmax)

        self.numgauss = QSpinBox()
        self.numgauss.setRange(10, 1000)
        self.numgauss.setValue(60)  # Default: 1.5 * 40 = 60
        self.numgauss.setToolTip("Number of Gauss-Laguerre points\nDefault: 1.5 × Rmax")
        left_form.addRow("No. of Gauss points:", self.numgauss)

        grid.addLayout(left_form, 0, 0)

        # Right column
        right_form = QFormLayout()

        self.rmaxgauss = QDoubleSpinBox()
        self.rmaxgauss.setRange(1, 500)
        self.rmaxgauss.setValue(40)  # Default: same as Rmax
        rmaxgauss_label = self.create_label_with_info(
            "Gauss Rmax (fm):",
            "Maximum radius for Gauss-Laguerre mesh (fm)\n"
            "Default: same as Rmax\n"
            "Should cover the entire interaction region"
        )
        right_form.addRow(rmaxgauss_label, self.rmaxgauss)

        self.method = QSpinBox()
        self.method.setRange(1, 2)
        self.method.setValue(1)
        self.method.valueChanged.connect(self.on_method_changed)
        method_label = self.create_label_with_info(
            "Method:",
            "Solution method for scattering problem:\n\n"
            "1 = Linear Equation Method\n"
            "    • Solves scattering directly from coupled equations\n"
            "    • Works with nonlocal potentials\n"
            "    • Generally more stable\n\n"
            "2 = Green's Function Method\n"
            "    • Solves bound states first, then uses Green's function\n"
            "    • Does NOT work with nonlocal potentials\n"
            "    • Can be faster for some problems\n\n"
            "Recommendation: Use Method 1 (default)"
        )
        right_form.addRow(method_label, self.method)

        self.cwftype = QSpinBox()
        self.cwftype.setRange(1, 2)
        self.cwftype.setValue(1)
        cwf_label = self.create_label_with_info(
            "Coulomb function type:",
            "Coulomb Wave Function calculation library:\n\n"
            "1 = COULCC (Fortran) ⭐ RECOMMENDED\n"
            "    • Fortran implementation using Steed's method\n"
            "    • Stable and well-tested\n"
            "    • MUCH MUCH faster performance\n"
            "    • Best choice for most calculations\n\n"
            "2 = C++ library (cwfcomplex)\n"
            "    • C++ implementation from adyo_v1_0\n"
            "    • Slower than COULCC\n"
            "    • Requires C++ library compilation\n\n"
            "⚡ Recommendation: Use 1 (default) - much faster!"
        )
        right_form.addRow(cwf_label, self.cwftype)

        self.thetah = QDoubleSpinBox()
        self.thetah.setRange(0.1, 10)
        self.thetah.setDecimals(2)
        self.thetah.setValue(1.0)
        thetah_label = self.create_label_with_info(
            "Scattering angle step (deg):",
            "Scattering angle step size (in degrees)\n"
            "Defines angular resolution for differential cross section σ(θ)\n"
            "Smaller step = finer angular grid, more data points\n"
            "Typical: 0.5° - 1.0°"
        )
        right_form.addRow(thetah_label, self.thetah)

        self.thetamax = QDoubleSpinBox()
        self.thetamax.setRange(0, 180)
        self.thetamax.setValue(180)
        thetamax_label = self.create_label_with_info(
            "Max scattering angle (deg):",
            "Maximum scattering angle (in degrees)\n"
            "Calculates differential cross section from 0° to θ_max\n"
            "180° = full angular range (forward to backward)\n"
            "90° = only forward hemisphere"
        )
        right_form.addRow(thetamax_label, self.thetamax)

        # Checkboxes in right column
        self.matgauss = QCheckBox()
        self.matgauss.setChecked(False)
        matgauss_label = self.create_label_with_info(
            "Use Gauss for matrix:",
            "Use Gauss-Laguerre quadrature for potential matrix elements\n\n"
            "What it does:\n"
            "  • True: Calculate full matrix <i|V|j> using Gauss quadrature\n"
            "          (accurate off-diagonal elements)\n"
            "  • False: Diagonal approximation <i|V|i> only\n"
            "          (faster but less accurate)\n\n"
            "Recommendation: True for accurate results\n"
            "Use False only for quick tests"
        )
        right_form.addRow(matgauss_label, self.matgauss)

        self.bgauss = QCheckBox()
        self.bgauss.setChecked(False)
        self.bgauss.stateChanged.connect(self.on_bgauss_changed)
        bgauss_label = self.create_label_with_info(
            "Use Gauss for B vector:",
            "Use Gauss-Laguerre quadrature for inhomogeneous term (B vector)\n\n"
            "What it does:\n"
            "  • True: Calculate B = ∫ φ(r) V(r) f_L(r) dr using Gauss quadrature\n"
            "          (accurate integration)\n"
            "  • False: Use Lagrange mesh points for integration\n"
            "          (faster but less accurate)\n\n"
            "NOTE: Nonlocal potentials REQUIRE this to be True!\n\n"
            "Recommendation: True for accurate results"
        )
        right_form.addRow(bgauss_label, self.bgauss)

        self.backrot = QCheckBox()
        self.backrot.setChecked(False)
        backrot_label = self.create_label_with_info(
            "Back-rotation:",
            "Back-rotation technique for complex scaling\n\n"
            "What it does:\n"
            "  • Rotates the coordinate back to real axis before calculating observables\n"
            "  • Improves numerical stability for cross sections\n"
            "  • Essential when using complex contour rotation (ctheta > 0)\n\n"
            "When to use:\n"
            "  • Always enable with contour rotation (ctheta > 0)\n"
            "  • Automatically enabled when reading external potential\n\n"
            "Recommendation: True (better stability)"
        )
        right_form.addRow(backrot_label, self.backrot)

        # Read external potential with file path display
        readinpot_layout = QHBoxLayout()
        self.readinpot = QCheckBox()
        self.readinpot.setChecked(False)
        self.readinpot.stateChanged.connect(self.on_readinpot_changed)
        readinpot_layout.addWidget(self.readinpot)

        self.pot_file_label = QLabel("(No file selected)")
        self.pot_file_label.setStyleSheet("color: gray; font-size: 10px;")
        readinpot_layout.addWidget(self.pot_file_label)
        readinpot_layout.addStretch()

        readinpot_label = self.create_label_with_info(
            "Read external potential:",
            "Read optical potential from external file\n\n"
            "What it does:\n"
            "  • Loads pre-calculated potential V(r) from data file\n"
            "  • File should contain radial mesh and potential values\n"
            "  • Useful for microscopic potentials from structure calculations\n\n"
            "When enabled:\n"
            "  • Opens file dialog to select potential file (.dat or .txt)\n"
            "  • Automatically enables back-rotation for stability\n"
            "  • File is copied to working directory before calculation\n\n"
            "File format: Two columns [radius(fm), V(r)]"
        )
        right_form.addRow(readinpot_label, readinpot_layout)

        grid.addLayout(right_form, 0, 1)

        layout.addLayout(grid)
        layout.addStretch()

        # Navigation buttons
        nav_layout = QHBoxLayout()
        nav_layout.addStretch()
        next_btn = QPushButton("Next: System Parameters →")
        next_btn.clicked.connect(lambda: self.tabs.setCurrentIndex(1))
        nav_layout.addWidget(next_btn)
        layout.addLayout(nav_layout)

        return widget

    def create_system_tab(self):
        """Create system parameters tab"""
        widget = QWidget()
        layout = QVBoxLayout(widget)

        # Title
        title = QLabel("System Parameters")
        title.setStyleSheet("font-size: 16px; font-weight: bold;")
        title.setAlignment(Qt.AlignCenter)
        layout.addWidget(title)

        # Two-column grid
        grid = QGridLayout()
        grid.setColumnStretch(0, 1)
        grid.setColumnStretch(1, 1)

        # Left column - Projectile
        left_form = QFormLayout()

        self.namep = QLineEdit("n")
        self.namep.setToolTip("Projectile name")
        left_form.addRow("Projectile Name:", self.namep)

        self.massp = QDoubleSpinBox()
        self.massp.setRange(0, 300)
        self.massp.setDecimals(1)
        self.massp.setValue(1)
        self.massp.setToolTip("Projectile mass (amu)")
        left_form.addRow("Projectile Mass (amu):", self.massp)

        self.zp = QDoubleSpinBox()
        self.zp.setRange(-10, 120)
        self.zp.setDecimals(1)
        self.zp.setValue(0)
        self.zp.setToolTip("Projectile charge")
        left_form.addRow("Projectile Charge:", self.zp)

        self.sp = QDoubleSpinBox()
        self.sp.setRange(0, 10)
        self.sp.setDecimals(1)
        self.sp.setValue(0.5)
        self.sp.setToolTip("Projectile spin")
        left_form.addRow("Projectile Spin:", self.sp)

        grid.addLayout(left_form, 0, 0)

        # Right column - Target
        right_form = QFormLayout()

        self.namet = QLineEdit("40Ca")
        self.namet.setToolTip("Target name")
        right_form.addRow("Target Name:", self.namet)

        self.masst = QDoubleSpinBox()
        self.masst.setRange(0, 300)
        self.masst.setDecimals(1)
        self.masst.setValue(40)
        self.masst.setToolTip("Target mass (amu)")
        right_form.addRow("Target Mass (amu):", self.masst)

        self.zt = QDoubleSpinBox()
        self.zt.setRange(-10, 120)
        self.zt.setDecimals(1)
        self.zt.setValue(20)
        self.zt.setToolTip("Target charge")
        right_form.addRow("Target Charge:", self.zt)

        self.elab = QDoubleSpinBox()
        self.elab.setRange(0, 1000)
        self.elab.setDecimals(1)
        self.elab.setValue(20)
        self.elab.setToolTip("Laboratory energy (MeV)")
        right_form.addRow("Lab Energy (MeV):", self.elab)

        self.jmin = QDoubleSpinBox()
        self.jmin.setRange(0, 200)
        self.jmin.setDecimals(1)
        self.jmin.setValue(0)
        self.jmin.setToolTip("Minimum total angular momentum")
        right_form.addRow("Min. Angular Mom. (J):", self.jmin)

        self.jmax = QDoubleSpinBox()
        self.jmax.setRange(0, 200)
        self.jmax.setDecimals(1)
        self.jmax.setValue(10)
        self.jmax.setToolTip("Maximum total angular momentum")
        right_form.addRow("Max. Angular Mom. (J):", self.jmax)

        grid.addLayout(right_form, 0, 1)

        layout.addLayout(grid)
        layout.addStretch()

        # Navigation buttons
        nav_layout = QHBoxLayout()
        prev_btn = QPushButton("← Previous: General")
        prev_btn.clicked.connect(lambda: self.tabs.setCurrentIndex(0))
        nav_layout.addWidget(prev_btn)
        nav_layout.addStretch()
        next_btn = QPushButton("Next: Potential →")
        next_btn.clicked.connect(lambda: self.tabs.setCurrentIndex(2))
        nav_layout.addWidget(next_btn)
        layout.addLayout(nav_layout)

        return widget

    def create_potential_tab(self):
        """Create optical potential parameters tab"""
        widget = QWidget()
        main_layout = QVBoxLayout(widget)
        main_layout.setContentsMargins(0, 0, 0, 0)

        # Create scroll area for potential parameters
        from PySide6.QtWidgets import QScrollArea
        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOff)

        # Container for all content
        container = QWidget()
        layout = QVBoxLayout(container)

        # Title
        title = QLabel("Optical Model Potential")
        title.setStyleSheet("font-size: 16px; font-weight: bold;")
        title.setAlignment(Qt.AlignCenter)
        layout.addWidget(title)

        # Volume Terms
        volume_group = QGroupBox("Volume Terms")
        volume_layout = QGridLayout()

        # Real Volume (left)
        real_vol = QLabel("Real Volume")
        real_vol.setStyleSheet("font-weight: bold;")
        volume_layout.addWidget(real_vol, 0, 0)

        real_vol_form = QFormLayout()
        self.vv = QDoubleSpinBox()
        self.vv.setRange(-500, 500)
        self.vv.setDecimals(3)
        self.vv.setValue(46.553)
        self.vv.setToolTip("Depth of real volume potential (MeV)")
        real_vol_form.addRow("Depth (MeV):", self.vv)

        self.rv = QDoubleSpinBox()
        self.rv.setRange(0, 10)
        self.rv.setDecimals(3)
        self.rv.setValue(1.185)
        self.rv.setToolTip("Radius parameter (fm)")
        real_vol_form.addRow("Radius (fm):", self.rv)

        self.av = QDoubleSpinBox()
        self.av.setRange(0, 5)
        self.av.setDecimals(3)
        self.av.setValue(0.672)
        self.av.setToolTip("Diffuseness (fm)")
        real_vol_form.addRow("Diffuseness (fm):", self.av)

        volume_layout.addLayout(real_vol_form, 1, 0)

        # Imaginary Volume (right)
        imag_vol = QLabel("Imaginary Volume")
        imag_vol.setStyleSheet("font-weight: bold;")
        volume_layout.addWidget(imag_vol, 0, 1)

        imag_vol_form = QFormLayout()
        self.wv = QDoubleSpinBox()
        self.wv.setRange(-500, 500)
        self.wv.setDecimals(3)
        self.wv.setValue(1.777)
        self.wv.setToolTip("Depth of imaginary volume potential (MeV)")
        imag_vol_form.addRow("Depth (MeV):", self.wv)

        self.rw = QDoubleSpinBox()
        self.rw.setRange(0, 10)
        self.rw.setDecimals(3)
        self.rw.setValue(1.185)
        self.rw.setToolTip("Radius parameter (fm)")
        imag_vol_form.addRow("Radius (fm):", self.rw)

        self.aw = QDoubleSpinBox()
        self.aw.setRange(0, 5)
        self.aw.setDecimals(3)
        self.aw.setValue(0.672)
        self.aw.setToolTip("Diffuseness (fm)")
        imag_vol_form.addRow("Diffuseness (fm):", self.aw)

        volume_layout.addLayout(imag_vol_form, 1, 1)
        volume_group.setLayout(volume_layout)
        layout.addWidget(volume_group)

        # Surface Terms
        surface_group = QGroupBox("Surface Terms")
        surface_layout = QGridLayout()

        # Real Surface (left)
        real_surf = QLabel("Real Surface")
        real_surf.setStyleSheet("font-weight: bold;")
        surface_layout.addWidget(real_surf, 0, 0)

        real_surf_form = QFormLayout()
        self.vs = QDoubleSpinBox()
        self.vs.setRange(-500, 500)
        self.vs.setDecimals(3)
        self.vs.setValue(0)
        self.vs.setToolTip("Depth of real surface potential (MeV)")
        real_surf_form.addRow("Depth (MeV):", self.vs)

        self.rvs = QDoubleSpinBox()
        self.rvs.setRange(0, 10)
        self.rvs.setDecimals(3)
        self.rvs.setValue(0)
        self.rvs.setToolTip("Radius parameter (fm)")
        real_surf_form.addRow("Radius (fm):", self.rvs)

        self.avs = QDoubleSpinBox()
        self.avs.setRange(0, 5)
        self.avs.setDecimals(3)
        self.avs.setValue(0)
        self.avs.setToolTip("Diffuseness (fm)")
        real_surf_form.addRow("Diffuseness (fm):", self.avs)

        surface_layout.addLayout(real_surf_form, 1, 0)

        # Imaginary Surface (right)
        imag_surf = QLabel("Imaginary Surface")
        imag_surf.setStyleSheet("font-weight: bold;")
        surface_layout.addWidget(imag_surf, 0, 1)

        imag_surf_form = QFormLayout()
        self.ws = QDoubleSpinBox()
        self.ws.setRange(-500, 500)
        self.ws.setDecimals(3)
        self.ws.setValue(7.182)
        self.ws.setToolTip("Depth of imaginary surface potential (MeV)")
        imag_surf_form.addRow("Depth (MeV):", self.ws)

        self.rws = QDoubleSpinBox()
        self.rws.setRange(0, 10)
        self.rws.setDecimals(3)
        self.rws.setValue(1.288)
        self.rws.setToolTip("Radius parameter (fm)")
        imag_surf_form.addRow("Radius (fm):", self.rws)

        self.aws = QDoubleSpinBox()
        self.aws.setRange(0, 5)
        self.aws.setDecimals(3)
        self.aws.setValue(0.538)
        self.aws.setToolTip("Diffuseness (fm)")
        imag_surf_form.addRow("Diffuseness (fm):", self.aws)

        surface_layout.addLayout(imag_surf_form, 1, 1)
        surface_group.setLayout(surface_layout)
        layout.addWidget(surface_group)

        # Spin-Orbit Terms
        so_group = QGroupBox("Spin-Orbit Terms")
        so_layout = QGridLayout()

        # Real SO (left)
        real_so = QLabel("Real Spin-Orbit")
        real_so.setStyleSheet("font-weight: bold;")
        so_layout.addWidget(real_so, 0, 0)

        real_so_form = QFormLayout()
        self.vsov = QDoubleSpinBox()
        self.vsov.setRange(-500, 500)
        self.vsov.setDecimals(3)
        self.vsov.setValue(5.343)
        self.vsov.setToolTip("Depth of real spin-orbit potential (MeV)")
        real_so_form.addRow("Depth (MeV):", self.vsov)

        self.rsov = QDoubleSpinBox()
        self.rsov.setRange(0, 10)
        self.rsov.setDecimals(3)
        self.rsov.setValue(0.996)
        self.rsov.setToolTip("Radius parameter (fm)")
        real_so_form.addRow("Radius (fm):", self.rsov)

        self.asov = QDoubleSpinBox()
        self.asov.setRange(0, 5)
        self.asov.setDecimals(3)
        self.asov.setValue(0.590)
        self.asov.setToolTip("Diffuseness (fm)")
        real_so_form.addRow("Diffuseness (fm):", self.asov)

        so_layout.addLayout(real_so_form, 1, 0)

        # Imaginary SO (right)
        imag_so = QLabel("Imaginary Spin-Orbit")
        imag_so.setStyleSheet("font-weight: bold;")
        so_layout.addWidget(imag_so, 0, 1)

        imag_so_form = QFormLayout()
        self.vsow = QDoubleSpinBox()
        self.vsow.setRange(-500, 500)
        self.vsow.setDecimals(3)
        self.vsow.setValue(-0.110)
        self.vsow.setToolTip("Depth of imaginary spin-orbit potential (MeV)")
        imag_so_form.addRow("Depth (MeV):", self.vsow)

        self.rsow = QDoubleSpinBox()
        self.rsow.setRange(0, 10)
        self.rsow.setDecimals(3)
        self.rsow.setValue(0.996)
        self.rsow.setToolTip("Radius parameter (fm)")
        imag_so_form.addRow("Radius (fm):", self.rsow)

        self.asow = QDoubleSpinBox()
        self.asow.setRange(0, 5)
        self.asow.setDecimals(3)
        self.asow.setValue(0.590)
        self.asow.setToolTip("Diffuseness (fm)")
        imag_so_form.addRow("Diffuseness (fm):", self.asow)

        so_layout.addLayout(imag_so_form, 1, 1)
        so_group.setLayout(so_layout)
        layout.addWidget(so_group)

        # Coulomb Term
        coulomb_group = QGroupBox("Coulomb Term & Radius Parameters")
        coulomb_form = QFormLayout()

        self.rc = QDoubleSpinBox()
        self.rc.setRange(0, 10)
        self.rc.setDecimals(3)
        self.rc.setValue(1.698)
        self.rc.setToolTip("Coulomb radius parameter (fm)")
        coulomb_form.addRow("Coulomb rc (fm):", self.rc)

        a1_label = self.create_label_with_info(
            "a1:",
            "Mass number of projectile for potential radius calculation\n\n"
            "What it does:\n"
            "  • Used to calculate R₀ = r₀ × (a1^(1/3) + a2^(1/3))\n"
            "  • R₀ is the potential radius parameter\n\n"
            "Default behavior (a1=0, a2=0):\n"
            "  • a1 = 0\n"
            "  • a2 = target mass (masst)\n"
            "  • R₀ = r₀ × masst^(1/3)\n\n"
            "When to set:\n"
            "  • For asymmetric systems (e.g., p + ²⁰⁸Pb)\n"
            "  • Leave at 0 for default behavior (recommended)\n\n"
            "Typical: 0 (use default)"
        )
        self.a1 = QDoubleSpinBox()
        self.a1.setRange(0, 300)
        self.a1.setDecimals(2)
        self.a1.setValue(0.0)
        self.a1.setToolTip("Mass number of projectile (0 = use default)")
        coulomb_form.addRow(a1_label, self.a1)

        a2_label = self.create_label_with_info(
            "a2:",
            "Mass number of target for potential radius calculation\n\n"
            "What it does:\n"
            "  • Used to calculate R₀ = r₀ × (a1^(1/3) + a2^(1/3))\n"
            "  • R₀ is the potential radius parameter\n\n"
            "Default behavior (a1=0, a2=0):\n"
            "  • a1 = 0\n"
            "  • a2 = target mass (masst)\n"
            "  • R₀ = r₀ × masst^(1/3)\n\n"
            "When to set:\n"
            "  • For asymmetric systems (e.g., p + ²⁰⁸Pb)\n"
            "  • Leave at 0 for default behavior (recommended)\n\n"
            "Typical: 0 (use default)"
        )
        self.a2 = QDoubleSpinBox()
        self.a2.setRange(0, 300)
        self.a2.setDecimals(2)
        self.a2.setValue(0.0)
        self.a2.setToolTip("Mass number of target (0 = use default)")
        coulomb_form.addRow(a2_label, self.a2)

        coulomb_group.setLayout(coulomb_form)
        layout.addWidget(coulomb_group)

        layout.addStretch()

        # Set container in scroll area
        scroll.setWidget(container)
        main_layout.addWidget(scroll)

        # Navigation buttons (outside scroll area)
        nav_layout = QHBoxLayout()
        prev_btn = QPushButton("← Previous: System")
        prev_btn.clicked.connect(lambda: self.tabs.setCurrentIndex(1))
        nav_layout.addWidget(prev_btn)
        nav_layout.addStretch()
        next_btn = QPushButton("Next: Nonlocal →")
        next_btn.clicked.connect(lambda: self.tabs.setCurrentIndex(3))
        nav_layout.addWidget(next_btn)
        main_layout.addLayout(nav_layout)

        return widget

    def create_nonlocal_tab(self):
        """Create nonlocal potential parameters tab"""
        widget = QWidget()
        layout = QVBoxLayout(widget)

        # Title
        title = QLabel("Nonlocal Potential Parameters")
        title.setStyleSheet("font-size: 16px; font-weight: bold;")
        title.setAlignment(Qt.AlignCenter)
        layout.addWidget(title)

        # Form
        form = QFormLayout()

        self.use_nonlocal = QCheckBox()
        self.use_nonlocal.setChecked(False)
        self.use_nonlocal.setToolTip("Enable nonlocal potential\n\n"
                                     "NOTE: Nonlocal potential requires Method 1 (Linear Equation)\n"
                                     "Method 2 (Green's Function) does NOT support nonlocal potentials")
        self.use_nonlocal.stateChanged.connect(self.on_nonlocal_changed)
        form.addRow("Use nonlocal potential:", self.use_nonlocal)

        self.nlbeta = QDoubleSpinBox()
        self.nlbeta.setRange(0, 10)
        self.nlbeta.setDecimals(2)
        self.nlbeta.setValue(0.0)
        self.nlbeta.setToolTip("Nonlocal range parameter (fm)")
        form.addRow("Nonlocal range (fm):", self.nlbeta)

        layout.addLayout(form)

        # Info text
        info = QLabel("The nonlocal form follows Perey and Buck (Nuclear Physics, 32, 353-380).\n"
                     "If enabled, beta value will be used as the nonlocal range parameter.")
        info.setStyleSheet("color: gray; font-size: 11px;")
        info.setWordWrap(True)
        layout.addWidget(info)

        layout.addStretch()

        # Navigation buttons
        nav_layout = QHBoxLayout()
        prev_btn = QPushButton("← Previous: Potential")
        prev_btn.clicked.connect(lambda: self.tabs.setCurrentIndex(2))
        nav_layout.addWidget(prev_btn)
        nav_layout.addStretch()
        layout.addLayout(nav_layout)

        return widget

    def on_rmax_changed(self, value):
        """Update dependent values when Rmax changes"""
        # Update rmaxgauss to match rmax
        self.rmaxgauss.setValue(value)
        # Update numgauss to 1.5 * rmax (rounded to integer)
        self.numgauss.setValue(int(1.5 * value))

    def on_method_changed(self, value):
        """Warn if Method 2 is selected with nonlocal potential"""
        if value == 2 and self.use_nonlocal.isChecked():
            from PySide6.QtWidgets import QMessageBox
            QMessageBox.warning(
                self,
                "Invalid Configuration",
                "Method 2 (Green's Function) does NOT work with nonlocal potentials!\n\n"
                "Please either:\n"
                "• Switch to Method 1 (Linear Equation), or\n"
                "• Disable nonlocal potential"
            )
            # Auto-switch back to Method 1
            self.method.setValue(1)

    def on_bgauss_changed(self, state):
        """Warn if bgauss is disabled with nonlocal potential"""
        if state == Qt.CheckState.Unchecked.value and self.use_nonlocal.isChecked():
            from PySide6.QtWidgets import QMessageBox
            QMessageBox.warning(
                self,
                "Invalid Configuration",
                "Nonlocal potentials REQUIRE 'Use Gauss for B vector' to be enabled!\n\n"
                "Automatically re-enabling it."
            )
            # Auto-enable bgauss
            self.bgauss.setChecked(True)

    def on_nonlocal_changed(self, state):
        """Handle nonlocal potential checkbox change"""
        if state == Qt.CheckState.Checked.value:
            # Check Method 2 conflict
            if self.method.value() == 2:
                from PySide6.QtWidgets import QMessageBox
                QMessageBox.warning(
                    self,
                    "Invalid Configuration",
                    "Nonlocal potentials do NOT work with Method 2 (Green's Function)!\n\n"
                    "Automatically switching to Method 1 (Linear Equation)."
                )
                # Auto-switch to Method 1
                self.method.setValue(1)

            # Auto-enable bgauss (required for nonlocal)
            if not self.bgauss.isChecked():
                self.bgauss.setChecked(True)

    def on_readinpot_changed(self, state):
        """Handle read external potential checkbox change"""
        if state == Qt.CheckState.Checked.value:
            # Open file dialog to select potential file
            file_path, _ = QFileDialog.getOpenFileName(
                self,
                "Select Potential File",
                "",
                "Data Files (*.dat *.txt);;All Files (*)"
            )

            if file_path:
                # File selected - store path and update label
                self.pot_file_path = file_path
                # Show only filename in label
                import os
                filename = os.path.basename(file_path)
                self.pot_file_label.setText(f"({filename})")
                self.pot_file_label.setStyleSheet("color: green; font-size: 10px;")

                # Auto-check back-rotation
                self.backrot.setChecked(True)
            else:
                # No file selected - uncheck the checkbox
                self.readinpot.setChecked(False)
                self.pot_file_path = None
                self.pot_file_label.setText("(No file selected)")
                self.pot_file_label.setStyleSheet("color: gray; font-size: 10px;")
        else:
            # Unchecked - clear file path
            self.pot_file_path = None
            self.pot_file_label.setText("(No file selected)")
            self.pot_file_label.setStyleSheet("color: gray; font-size: 10px;")

    def generate_input(self):
        """Generate COLOSS input file content"""
        # Helper function for boolean conversion
        def bool_to_fortran(val):
            return 't' if val else 'f'

        # Note: method field doesn't exist in original code, so we'll skip it
        input_text = f"""&general
    nr={self.nr.value()}  alpha={self.alpha.value()} Rmax={self.rmax.value()} ctheta={self.ctheta.value()}
    matgauss={bool_to_fortran(self.matgauss.isChecked())} bgauss={bool_to_fortran(self.bgauss.isChecked())}
    numgauss={self.numgauss.value()} rmaxgauss={self.rmaxgauss.value()}
    thetah={self.thetah.value()} thetamax={self.thetamax.value()}
    cwftype={self.cwftype.value()} backrot={bool_to_fortran(self.backrot.isChecked())}
    readinpot={bool_to_fortran(self.readinpot.isChecked())} /

&system
    zp={self.zp.value()}    massp={self.massp.value()}   namep='{self.namep.text()}'
    zt={self.zt.value()}   masst={self.masst.value()}  namet='{self.namet.text()}'
    jmin={self.jmin.value()} jmax={self.jmax.value()}  elab={self.elab.value()} sp={self.sp.value()}/

&pot
    a1={self.a1.value():.2f}  a2={self.a2.value():.2f}
    vv={self.vv.value():.3f}   rv={self.rv.value():.3f}    av={self.av.value():.3f}
    wv={self.wv.value():.3f}    rw={self.rw.value():.3f}    aw={self.aw.value():.3f}
    vs={self.vs.value():.3f}        rvs={self.rvs.value():.3f}       avs={self.avs.value():.3f}
    ws={self.ws.value():.3f}    rws={self.rws.value():.3f}   aws={self.aws.value():.3f}
    vsov={self.vsov.value():.3f}  rsov={self.rsov.value():.3f}  asov={self.asov.value():.3f}
    vsow={self.vsow.value():.3f} rsow={self.rsow.value():.3f}  asow={self.asow.value():.3f}
    rc={self.rc.value():.3f} /

&nonlocalpot
    nonlocal={bool_to_fortran(self.use_nonlocal.isChecked())} nlbeta={self.nlbeta.value()} /
"""
        self.input_changed.emit()
        return input_text

    def get_input_text(self):
        """Get the generated input text"""
        return self.generate_input()

    def get_pot_file_path(self):
        """Get the path to the external potential file"""
        return self.pot_file_path
