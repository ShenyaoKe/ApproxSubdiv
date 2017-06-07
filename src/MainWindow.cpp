#include "MainWindow.h"

MainWindow::MainWindow(QWidget* parent)
	: QMainWindow(parent)
	, mUpdatePending(false), mAnimating(false)
	, mOglViewer(new OglViewer)
{

	mUI.setupUi(this);
	mUI.ogl_layout->addWidget(mOglViewer);
	//setWindowTitle(tr("OpenGL Qt Template"));

	mOglViewer->setFocusPolicy(Qt::StrongFocus);
}

MainWindow::~MainWindow()
{
	delete mOglViewer;
}

void MainWindow::triggerAboutWindow()
{
	mAbout = new QDialog(0,0);
	Ui::about_dialog about_ui;
	about_ui.setupUi(mAbout);
	mAbout->show();
}

void MainWindow::closeEvent(QCloseEvent *e)
{
	Q_UNUSED(e);

	foreach(QWidget *widget, QApplication::topLevelWidgets())
	{
		widget->close();
	}
}
