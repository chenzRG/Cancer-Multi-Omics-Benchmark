import torch

class Earlystopping:
    def __init__(self, number=50, path='../ssd/checkpoint.pt'):
        self.number = number
        self.path = path
        self.counter = 0
        self.epoch_count = 0
        self.best_epoch_num = 1
        self.max_acc = None
        self.stop_now = False

    def __call__(self, model, acc):
        if self.max_acc is None:
            self.epoch_count += 1
            self.max_acc = acc
            self.save_model(model)
        elif acc > self.max_acc:
            self.epoch_count += 1
            self.best_epoch_num = self.epoch_count
            print('New maximum accuracy: {:.4f}% -> {:.4f}%\n'.format(self.max_acc, acc))
            self.max_acc = acc
            self.save_model(model)
            self.counter = 0
        else:
            self.epoch_count += 1
            print('Current maximum accuracy: {:.4f}%\n'
                  'Early stopping counter: {:d}/{:d}\n'.format(self.max_acc, self.counter, self.number))
            if self.counter >= self.number:
                self.stop_now = True
            self.counter += 1

    def save_model(self, model):
        torch.save(model.state_dict(), self.path)
