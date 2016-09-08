#ifndef DATASELECTDIALOG_H
#define DATASELECTDIALOG_H

//forward declarations
struct InputParameters;

#include <QDialog>
#include <QFile>
#include <QString>
#include <QtWidgets>

#include <map>

namespace Ui {
class DataSelectDialog;
}

class DataSelectDialog : public QDialog
{
    Q_OBJECT

public:
    explicit DataSelectDialog(InputParameters& params, QWidget *parent = 0);
    ~DataSelectDialog();

signals:
    void dataSelected();

protected:
    void closeEvent(QCloseEvent* event);

private slots:
    void on_computeButton_clicked();
    void on_openFileButton_clicked();

//TESTING ONLY
//    void on_shortcutButton1_clicked();
//    void on_shortcutButton2_clicked();

private:
    Ui::DataSelectDialog *ui;

    InputParameters& params;

    void invalid_file(const QString& message);
    void detect_file_type();

    bool data_selected;     //false until user clicks "compute" button

};

#endif // DATASELECTDIALOG_H
