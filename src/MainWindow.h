#pragma once
#ifndef QTOGL_WINDOW_H
#define QTOGL_WINDOW_H

#include <QMainWindow>
#include <QDialog>

#include "OglViewer.h"

#include "ui_MainWindow.h"
#include "ui_About.h"


class MainWindow : public QMainWindow
{
	Q_OBJECT

public:
	MainWindow(QWidget* parent = 0);
	~MainWindow();

public:
	void triggerAboutWindow();
	//void aboutwindow();
protected:
	void closeEvent(QCloseEvent* e);
private:
	bool mUpdatePending;
	bool mAnimating;

	OglViewer* mOglViewer;
	QDialog* mAbout;
	Ui::MainWindowClass mUI;
};

#endif // QTOGL_WINDOW_H
